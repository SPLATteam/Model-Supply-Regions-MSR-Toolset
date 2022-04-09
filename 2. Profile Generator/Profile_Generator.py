import pandas as pd
from netCDF4 import Dataset
import numpy as np
import itertools
from scipy import spatial
from scipy.interpolate import interp1d
from rasterstats import zonal_stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import pvlib
import re
import math
import geopandas as gpd
import time
import traceback

def ComputeHourlyCF_SolarPV (GSA_GHI_MSR_Mean, temp_2m, GHI): #returns hourly CF for the location whose hourly inputs are given as inputs
    # Inputs
        # GHI= GHI values (8760 element numpy array)
        # temp_2m= Temperature values at 2 meter (8760 element numpy array)
        # GSA_GHI_MSR_Mean= Global Solar Atlas daily mean available at 1km resolution averaged across the MSR
    # output:
        # CF_solar= Solar CFs (8760 element numpy array)
        # GHI_corrected_wh= Corrected GHI (8760 element numpy array)
        # BiasCorrection_GHI_Adder_Wh
        # ERA5_GHI_OriginalAnnualYield (for any later use e.g. debugging etc)

    GHI=np.ma.filled(GHI)
    temp_2m=np.ma.filled(temp_2m)
    GHI[GHI < 0] = 0
    GHI=GHI/3600000 #convert ERA GHI given in jouls/m2 to kwh/m2 (the units of GSA GHI dataset)
    ERA5_GHI_OriginalAnnualYield=np.sum(GHI)
    GHIAnnualBias_KWh=GSA_GHI_MSR_Mean*(8760/24)-ERA5_GHI_OriginalAnnualYield #GSA GHI dataset is mean per day yield
    GHI_ifallCorrected = GHI + GHIAnnualBias_KWh / 8760
    BiasCorrection_GHI_Adder_KWh=GHIAnnualBias_KWh/len(GHI[(GHI!=0)&(GHI_ifallCorrected<1)&(GHI_ifallCorrected>0)])
    GHI[(GHI!=0)&(GHI_ifallCorrected<1)&(GHI_ifallCorrected>0)]=GHI[(GHI!=0)&(GHI_ifallCorrected<1)&(GHI_ifallCorrected>0)]+BiasCorrection_GHI_Adder_KWh
    GHI[GHI<0]=0
    GHI_corrected_wh = GHI*1000
    BiasCorrection_GHI_Adder_Wh=BiasCorrection_GHI_Adder_KWh * 1000


    Gstc = 1000 #W/M^2
    Tmod_stc = 25 +273 #K
    k = [-0.017162, -0.040289, -0.004681, 0.000148, 0.000169, 0.000005]
    cT = 0.035 #per kelvin
    GHI_Wh=GHI_corrected_wh #Yunshu code uses W/m2
    G_norm = GHI_Wh/Gstc
    T_norm = cT*GHI_Wh + temp_2m-Tmod_stc
    log_G_norm = np.zeros(len(GHI_Wh))
    for n in range(len(GHI_Wh)):
        if GHI_Wh[n] == 0:
            log_G_norm[n] = 0
        else:
            log_G_norm[n] = np.log(G_norm[n])
    n_rel = 1 + k[0]*log_G_norm + k[1]*(log_G_norm)**2 + T_norm*(k[2] + k[3]*log_G_norm + k[4]*(log_G_norm)**2) + k[5]*T_norm**2
    CF0 = n_rel * G_norm
    CF_solar = np.nan_to_num(CF0)
    CF_solar[CF_solar<0]=0
    return CF_solar, GHI_corrected_wh, BiasCorrection_GHI_Adder_Wh, ERA5_GHI_OriginalAnnualYield

def ComputeHourlyCF_Wind(BiasCorrEffectiveWindSpeeds, u100, u10, v100, v10, WindTurbineHeight_meters, temp_2m, elevation, pd_WindSpeed_to_Power):# returns wind hourly CF
    #inputs
        #BiasCorrEffectiveWindSpeeds= 100 meter ERA5 Effective wind speeds bias corrected with Global Wind Atlas (8760 element numpy array)
        #u100, v100, u10, v10= ERA5 wind speed components needed to get wind speeds and resulting CFs at any desired height (8760 element numpy array)
        #WindTurbineHeight_meters= Desired height to get ERA5 resource profile (biascorrected at 100m) and CFs
        #temp_2m= Temperature values at 2 meter (8760 element numpy array)
        #elevation= MSR elevation got from ERA5, needed to get location's air density that would adjust wind speeds as proxy to impact of air density on wind CF
    #Ouput
        #CF_wind= Height and air density adjusted wind CFs (8760 element numpy array)
        #map_windspeed_hubh= effective wind speeds at desired height (it is not air density adjusted)
        #ERA5_AnnualMean100m: for later use or debugging purpose

    h_0 = 100.0 # Reference extrapolation height
    h_hub = WindTurbineHeight_meters

    #calculate pressure using IRENA-LBNL report, units mentioned there, temperature from ERA5 is in kelvins with no zeros
    rho_0 = 1.225  # kg/m^3
    R = 287.058  # J/(kg.K)
    pressure_ts = 101325 / np.exp(elevation * 9.80665 / (R * temp_2m))
    LocationDensity_ts = pressure_ts / (R * temp_2m)

    #calculate adjust wind speeds as per ACEC methodology
    map_windspeed_10m = (u10 ** 2 + v10 ** 2) ** (1 / 2)
    map_windspeed_100m = (u100 ** 2 + v100 ** 2) ** (1 / 2)
    ERA5_AnnualMean100m = np.mean(map_windspeed_100m)
    map_windspeed_100mCorrected = BiasCorrEffectiveWindSpeeds
    map_windspeed_100m_UnCorrected = map_windspeed_100m #only for diagnosis

    alpha = np.log(np.nan_to_num(map_windspeed_100m/map_windspeed_10m, posinf=0, neginf=0)) / np.log(100 / 10)#for alpha dont use corrected speeds
    if h_hub==h_0:
        map_windspeed_hubh=map_windspeed_100mCorrected
        map_windspeed_hubh_Uncorrected = map_windspeed_100m_UnCorrected #only for diagnosis
    else:
        map_windspeed_hubh = map_windspeed_100mCorrected * ((h_hub / h_0)** alpha)
        map_windspeed_hubh_Uncorrected = map_windspeed_100m_UnCorrected * ((h_hub / h_0) ** alpha) #only for diagnosis

    map_windspeed_hubh_adjusted = map_windspeed_hubh * (rho_0 / LocationDensity_ts) ** (1 / 3) #uplifts or downlifts the annual mean a bit from GWA


    win_avg_speed = map_windspeed_hubh.mean()
    win_class_bounds = [7.5, 8.5]
    if win_avg_speed <= win_class_bounds[0]:
        IEC_Class = 'IEC Class 3'
    elif win_class_bounds[0] <= win_avg_speed < win_class_bounds[1]:
        IEC_Class = 'IEC Class 2'
    elif win_avg_speed > win_class_bounds[1]:
        IEC_Class = 'IEC Class 1'
    else:
        raise Exception('unidentified speed detected')


    f = interp1d(pd_WindSpeed_to_Power.iloc[:, 0],pd_WindSpeed_to_Power[IEC_Class], fill_value = "extrapolate")
    # plt.figure(1)
    # plt.plot(np.arange(0, 25, 0.1),f(np.arange(0, 25, 0.1)))
    CF_wind = f(np.ma.filled(map_windspeed_hubh_adjusted[:]))

    return CF_wind, map_windspeed_hubh, map_windspeed_hubh_Uncorrected, ERA5_AnnualMean100m

def develop_allMSR_8760BiasCorrEffectiveWindSpeeds(pd_allMSR_Hourly100mEffectiveWindSpeeds, pd_Wind_GWA_MSR_Mean):
    #Inputs
        #pd_allMSR_Hourly100mEffectiveWindSpeeds: All MSR 8760 hourly timeseries, pandas dataframe, Columns are Hour numbers i.e. H1, H2...
        #pd_Wind_GWA_MSR_Mean=All MSR Global Wind atlas annual mean wind speed at 1km resolution averaged across each MSR
    #Outputs
        #pd_allMSR_TimeSeriesProfilesCorrected

    def f(x, a, b):
        return a * x + b
    if len(pd_allMSR_Hourly100mEffectiveWindSpeeds)==1:
        return pd_allMSR_Hourly100mEffectiveWindSpeeds*float(pd_Wind_GWA_MSR_Mean/pd_allMSR_Hourly100mEffectiveWindSpeeds.mean(axis=1))
    else:

        # custom steps if needed
        pd_allMSR_TimeSeriesProfiles = pd_allMSR_Hourly100mEffectiveWindSpeeds
        pd_allMSR_GWAmean = pd_Wind_GWA_MSR_Mean

        # Phase-1 Compute x, y, fit ax+b for each rank, extrapolate each rank

        # compute the dependent variable x i.e. annual ERA5 means
        pd_allMSR_ERA5mean = pd_allMSR_TimeSeriesProfiles.mean(axis=1)

        ranks = [i + 1 for i in range(8760)]  # 1,2...8760
        # sort time series in ascending order. Reformat the column names as rank numbers [1,2...8760]. Rank 1 goes to the smallest wind speed. Data column under each rank will be indpendent variable y
        pd_allMSR_TimeSortedProfiles = pd_allMSR_TimeSeriesProfiles.transpose().apply(np.sort,
                                                                                        axis=0).transpose().set_axis(ranks,
                                                                                                                     axis=1)

        pd_allMSR_TimeSortedExtraplolatedProfiles = pd.DataFrame().reindex_like(
            pd_allMSR_TimeSortedProfiles)  # copy fromat of the dataframe

        # activate when diagnosis needed
        # pd_LinearFitParam_ab=pd.DataFrame()

        #plt.figure(1)
        for r in range(1, len(pd_allMSR_TimeSeriesProfiles.iloc[0,
                              :8760]) + 1):  # rank wise loop, get rank wise fitted equations, rank wise extrapolate (bias-correction) Time Sorted data values
            print("rank %s: curve fitted, value extrapolated"%r)
            popt, pcov = curve_fit(f, pd_allMSR_ERA5mean, pd_allMSR_TimeSortedProfiles.iloc[:,
                                                           r - 1])  # get a and b by fitting the ax+b over all MSR datapoints for a given rank 'r'. Column calling here is by numerical order i.e. (0,8760) so needed to subtract 1

            # activate when diagnosis needed
            # pd_LinearFitParam_ab = pd_LinearFitParam_ab.append(pd.DataFrame({'rank': [r], 'a': [popt[0]], 'b': [popt[1]]}))
            #plt.plot(pd_allMSR_ERA5mean, f(pd_allMSR_ERA5mean, popt[0], popt[1]), label='data')

            # extrapolated/Biascorrect
            pd_allMSR_TimeSortedExtraplolatedProfiles.iloc[:, r - 1] = f(pd_allMSR_GWAmean, popt[0], popt[1])


        # activate when diagnosis needed
        # plt.figure(1)
        # plt.plot(pd_allMSR_TimeSortedExtraplolatedProfiles.iloc[0, :])
        # plt.figure(2)
        # plt.plot(pd_allMSR_TimeSortedProfiles.iloc[0, :])

        # Phase 2 Mapping
        # Map the original time series profiles
        pd_allMSR_TimeSeriesProfileRanks = pd_allMSR_TimeSeriesProfiles.rank(ascending=True, method='first',
                                                                               axis=1)  # Get ranks of each TimeSeries data value, smallest gets 1, highest gets 8760

        # reverse map the corrected time sorted profiles to time series profiles
        pd_allMSR_TimeSeriesProfilesCorrected = pd.DataFrame().reindex_like(
            pd_allMSR_TimeSeriesProfiles)  # copy dataset format
        for z in range(0, len(pd_allMSR_TimeSeriesProfiles)):  # MSR wise loop
            pd_allMSR_TimeSeriesProfilesCorrected.iloc[z, :] = pd_allMSR_TimeSeriesProfileRanks.iloc[z, :].map(
                pd_allMSR_TimeSortedExtraplolatedProfiles.iloc[z, :])

    return pd_allMSR_TimeSeriesProfilesCorrected

def CreateLocalTimeProfile(pd_UTC, pd_CountryUTC_offsets, country_withspaces):
    # Create Local Time Profile
    UTCHourTags = ['H%s' % i for i in range(1, 8761)]
    offset = pd_CountryUTC_offsets[pd_CountryUTC_offsets.Country == country_withspaces].Hours.iloc[0]
    LocalTimeHourTags = ['H%s' % (i % 8760 + 1) for i in range(offset, 8760 + offset)]
    ColRenameDictionary = dict(zip(UTCHourTags, LocalTimeHourTags))
    print("%s UTC Offset:%s" % (country_withspaces, offset))
    # print(ColRenameDictionary)

    pd_LocalTime = pd_UTC.rename(columns=ColRenameDictionary)
    pd_LocalTime = pd_LocalTime.set_index(pd_LocalTime.MSR_ID[:])
    pd_LocalTime = pd_LocalTime.drop(['MSR_ID'], axis=1)

    cols = pd_LocalTime.columns.to_list()
    if offset > 0:
        NewColOrder = cols[:-8760] + cols[-offset:] + cols[-8760:-offset]
        pd_LocalTime = pd_LocalTime[NewColOrder]
    if offset < 0:
        NewColOrder = cols[:-8760] + cols[-8760 - offset:] + cols[-8760:-8760 - offset]
        pd_LocalTime = pd_LocalTime[NewColOrder]

    return pd_LocalTime

#Start of the main program

#Read control input file
ControlPathsAndNames=pd.read_excel('ControlFile_ProfileGenerator.xlsx', sheet_name="PathsAndNames", index_col=0)
ControlConfigurations=pd.read_excel('ControlFile_ProfileGenerator.xlsx', sheet_name="Configurations", index_col=0)
pd_CountryUTC_offsets=pd.read_excel('ControlFile_ProfileGenerator.xlsx', sheet_name="CountryUTC_Offset_InUse")

#load paths
np_ERA5Data=Dataset(ControlPathsAndNames.loc["ERA5DataFilePath"][0])
pd_WindSpeed_to_Power = pd.read_csv(ControlPathsAndNames.loc["ThreeIEC_TurbinePowerCurves"][0])
AllCountries=pd.read_csv(ControlPathsAndNames.loc["FileAddress_CountryNamesList"][0],names=["Ct"])
Input_MSR_Folder=ControlPathsAndNames.loc["Input_MSR_Folder"][0]
OutputFolder_UTCProfiles=ControlPathsAndNames.loc["OutputFolder"][0]+"\\Results_UTC_Profiles"
OutputFolder_LocalTime=ControlPathsAndNames.loc["OutputFolder"][0]+"\\Results_LocalTimeProfiles"
ResourceRasterCarryingSubFolderName=ControlPathsAndNames.loc["ResourceRasterCarryingSubFolderName"][0]
MSR_DataCarryingSubFolderName=ControlPathsAndNames.loc["MSR_DataCarryingSubFolderName"][0]
SolarPVNameConvention=ControlPathsAndNames.loc["SolarPVNameConvention"][0]
WindNameConvention=ControlPathsAndNames.loc["WindNameConvention"][0]

#load configurations
RE_TechnologyList=[] # naming as per three technology names for which MSR creator code creates MSRs
ResourceRasterNameList=[]

if ControlConfigurations.loc["Run code for SolarPV"][0]==1:
    RE_TechnologyList.append(SolarPVNameConvention)
    ResourceRasterNameList.append(ControlPathsAndNames.loc["SolarPV_ResourceRasterName"][0])
if ControlConfigurations.loc["Run code for Wind"][0]==1:
    RE_TechnologyList.append(WindNameConvention)
    ResourceRasterNameList.append(ControlPathsAndNames.loc["Wind_ResourceRasterName"][0])
    WindTurbineHeight_meters = ControlConfigurations.loc["For wind, give wind turbine height (meters)"][0]

flag_Diagnosis= ControlConfigurations.loc["Produce hourly resource profiles for diagnostics"][0]

#code run time flags
flag_RunBiasCorrCode = 1  # Not a user option
pd_LogFile=pd.DataFrame()
DateTimeStamp = time.localtime()
DateTimeStamp = "%s%s%s%s%s%s" % (DateTimeStamp.tm_year, DateTimeStamp.tm_mon, DateTimeStamp.tm_mday, DateTimeStamp.tm_hour, DateTimeStamp.tm_min,DateTimeStamp.tm_sec)

if not os.path.isdir(OutputFolder_UTCProfiles):
    os.makedirs(OutputFolder_UTCProfiles)
if not os.path.isdir(OutputFolder_LocalTime):
    os.makedirs(OutputFolder_LocalTime)

for CountryCounter in range(0,len(AllCountries)):#country wise loop
    country_withspaces=AllCountries.Ct[CountryCounter]
    country = AllCountries.Ct[CountryCounter].replace(" ", "")
    MSR_CountryFolder=Input_MSR_Folder+r"\%s"%country

    OutputCountryFolder_UTC=OutputFolder_UTCProfiles + r"\%s" % country
    if not os.path.isdir(OutputCountryFolder_UTC):
        os.makedirs(OutputCountryFolder_UTC)

    OutputCountryFolder_LocalTime = OutputFolder_LocalTime + r"\%s" % country
    if not os.path.isdir(OutputCountryFolder_LocalTime):
        os.makedirs(OutputCountryFolder_LocalTime)

    for TechCounter in range (0,len(RE_TechnologyList)):
        RE=RE_TechnologyList[TechCounter]
        ResourceRasterName=ResourceRasterNameList[TechCounter]


        #read RE MSR information
        try:
            gpd_MSR_Attributes = gpd.read_file(MSR_CountryFolder + '\\' + MSR_DataCarryingSubFolderName + '\\' + RE + "_FinalMSRs.shp")
            # get stats of resource value across each MSR as Json dictionary (dc)
            dc_ResourceStatsAcrossMSR = zonal_stats(
                                    MSR_CountryFolder + '\\' + MSR_DataCarryingSubFolderName + '\\' + RE + "_FinalMSRs.shp",
                                    MSR_CountryFolder + '\\' + ResourceRasterCarryingSubFolderName + '\\' + ResourceRasterName + "_projected.tif",
                                    stats="count min mean max median sum")
            gpd_MSR_Attributes['Longitude'] =  gpd_MSR_Attributes.to_crs('EPSG:4326').centroid.x
            gpd_MSR_Attributes['Latitude'] =  gpd_MSR_Attributes.to_crs('EPSG:4326').centroid.y
            gpd_MSR_Attributes = gpd_MSR_Attributes.drop(['geometry'], axis=1)


            #find nearest neighbors of MSR coordinates in ERA5 data set
            ERA5Locations=list(itertools.product(np_ERA5Data.variables["latitude"][:],np_ERA5Data.variables["longitude"][:]))
            MSR_CentroidLocations=list(zip(gpd_MSR_Attributes.Latitude,gpd_MSR_Attributes.Longitude))
            ERA5Grid = spatial.cKDTree(ERA5Locations)
            dist, indexes = ERA5Grid.query(MSR_CentroidLocations)
            ERA5Locations_near_MSR_CentroidLocations= np.array(ERA5Locations)[indexes]
            gpd_MSR_Attributes["ERA5Latitude"]=ERA5Locations_near_MSR_CentroidLocations[:,0]
            gpd_MSR_Attributes["ERA5Longitude"]=ERA5Locations_near_MSR_CentroidLocations[:,1]
            np_lon = np_ERA5Data.variables["longitude"][:]
            np_lat = np_ERA5Data.variables["latitude"][:]


            #Declare variables to run subprograms that get capacity factor profiles from resource values for each MSR
            np_allMSR_HourlyCF=np.zeros((len(gpd_MSR_Attributes),8760))

            np_allMSR_HourlyWindSpeeds_corrected_ResultantVector = np.zeros((len(gpd_MSR_Attributes), 8760)) #m/s
            np_allMSR_HourlyWindSpeeds_uncorrected_ResultantVector = np.zeros((len(gpd_MSR_Attributes), 8760))  # m/s
            np_allMSR_HourlyGHI = np.zeros((len(gpd_MSR_Attributes), 8760)) #Wh
            np_allMSR_HourlyGHI_corrected_Wh=np.zeros((len(gpd_MSR_Attributes), 8760))

                #initialize datasets for wind bias correction
            np_allMSR_ERA5_AnnualMean=np.zeros(len(gpd_MSR_Attributes))
            np_allMSR_Wind_GWA_MSR_Mean = np.zeros(len(gpd_MSR_Attributes))

                #initialize datasets for solar bias correction
            np_allMSR_GHI_ERA5OriginalAnnualYield = np.zeros(len(gpd_MSR_Attributes))
            np_allMSR_GHI_GSA_MSR_Mean = np.zeros(len(gpd_MSR_Attributes))
            np_allMSR_BiasCorrection_GHI_Adder_Wh = np.zeros(len(gpd_MSR_Attributes))

            if flag_RunBiasCorrCode and RE==WindNameConvention:
                #get effective windspeeds, 8760 Bias Corrected
                #make all MSR timeseries profiles
                pd_allMSR_Hourly100mEffectiveWindSpeeds=pd.DataFrame()
                for MSR_Counter in range(0, len(gpd_MSR_Attributes)): #MSR wise loop
                    print("appending effective ERA5 speed dataset MSR %s"%MSR_Counter)
                    lat = ERA5Locations_near_MSR_CentroidLocations[MSR_Counter, 0]
                    lon = ERA5Locations_near_MSR_CentroidLocations[MSR_Counter, 1]
                    index_lat = np.where(np_lat == lat)[0][0]
                    index_lon = np.where(np_lon == lon)[0][0]
                    u=np_ERA5Data.variables["u100"][:,index_lat, index_lon]
                    v=np_ERA5Data.variables["v100"][:,index_lat, index_lon]
                    pd_allMSR_Hourly100mEffectiveWindSpeeds=pd_allMSR_Hourly100mEffectiveWindSpeeds.append (pd.DataFrame((u**2+v**2)**(1/2)).transpose(), ignore_index=True)
                pd_Wind_GWA_MSR_Mean = pd.DataFrame(dc_ResourceStatsAcrossMSR)['mean'].fillna(pd.DataFrame(dc_ResourceStatsAcrossMSR)['mean'].mean())#double checked, as long the MSR ids are numbered 0,1,2..., this command puts right mean to right table cell
                #get corrected 100m wind speeds, rows=MSRs, columns=hours
                pd_allMSR_Hourly100m8760BiasCorrEffectiveWindSpeeds=develop_allMSR_8760BiasCorrEffectiveWindSpeeds(pd_allMSR_Hourly100mEffectiveWindSpeeds, pd_Wind_GWA_MSR_Mean)


            for MSR_Counter in range(0,len(gpd_MSR_Attributes)): #MSR wise loop
                lat=ERA5Locations_near_MSR_CentroidLocations[MSR_Counter,0]
                lon=ERA5Locations_near_MSR_CentroidLocations[MSR_Counter,1]
                index_lat = np.where(np_lat == lat)[0][0]
                index_lon = np.where(np_lon == lon)[0][0]
                np_Hourly2meterTemperature = np_ERA5Data.variables["t2m"][:, index_lat, index_lon]

                if RE==SolarPVNameConvention:
                    np_HourlyGHI = np_ERA5Data.variables["ssrd"][:, index_lat, index_lon]
                    GSA_GHI_MSR_Mean=dc_ResourceStatsAcrossMSR[MSR_Counter]['mean']
                    if not GSA_GHI_MSR_Mean:
                        GSA_GHI_MSR_Mean=np.sum(np_HourlyGHI)*(24/8760)

                    np_HourlyCF, GHI_corrected_wh, BiasCorrection_GHI_Adder_Wh, ERA5_GHI_OriginalAnnualYield=ComputeHourlyCF_SolarPV(GSA_GHI_MSR_Mean,np_Hourly2meterTemperature,np_HourlyGHI)

                    np_HourlyCF=np.round(np_HourlyCF,decimals=3)
                    np_HourlyGHI = np.round(np_HourlyGHI/3600, decimals=3) #display in Watt Hours because the GHI biascorrect adder is also in Wh
                    np_HourlyGHI_corrected_Wh = np.round(GHI_corrected_wh, decimals=3) # In Watt Hours because the GHI biascorrect adder is also in Wh

                    np_allMSR_HourlyCF[MSR_Counter]=np_HourlyCF
                    np_allMSR_HourlyGHI[MSR_Counter] = np_HourlyGHI
                    np_allMSR_HourlyGHI_corrected_Wh[MSR_Counter]=np_HourlyGHI_corrected_Wh
                    np_allMSR_GHI_ERA5OriginalAnnualYield[MSR_Counter] = ERA5_GHI_OriginalAnnualYield
                    np_allMSR_GHI_GSA_MSR_Mean[MSR_Counter]= GSA_GHI_MSR_Mean*(8760/24) #converting per day yield of Solar GIS to annual
                    np_allMSR_BiasCorrection_GHI_Adder_Wh[MSR_Counter]=BiasCorrection_GHI_Adder_Wh

                if RE==WindNameConvention:
                    np_HourlyU100mEastward=np_ERA5Data.variables["u100"][:,index_lat, index_lon]
                    np_HourlyV100mNorthward = np_ERA5Data.variables["v100"][:, index_lat, index_lon]
                    np_HourlyU10mEastward = np_ERA5Data.variables["u10"][:, index_lat, index_lon]
                    np_HourlyV10mNorthward = np_ERA5Data.variables["v10"][:, index_lat, index_lon]
                    Wind_GWA_MSR_Mean=dc_ResourceStatsAcrossMSR[MSR_Counter]['mean']
                    BiasCorrEffectiveWindSpeeds= pd_allMSR_Hourly100m8760BiasCorrEffectiveWindSpeeds.iloc[MSR_Counter, :]
                    if not Wind_GWA_MSR_Mean:
                        Wind_GWA_MSR_Mean=pd.DataFrame(dc_ResourceStatsAcrossMSR)['mean'].mean()
                    elevation=np_ERA5Data.variables["z"][0, index_lat, index_lon]/9.80665 #see orography variable description in https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview

                    np_HourlyCF, np_HourlyWindSpeeds_ResultantVector, np_HourlyWindSpeeds_UnCorrected_ResultantVector, ERA5_AnnualMean=ComputeHourlyCF_Wind(BiasCorrEffectiveWindSpeeds, np_HourlyU100mEastward, np_HourlyU10mEastward, np_HourlyV100mNorthward, np_HourlyV10mNorthward, WindTurbineHeight_meters, np_Hourly2meterTemperature, elevation, pd_WindSpeed_to_Power)


                    np_HourlyCF = np.round(np_HourlyCF, decimals=3)
                    np_HourlyWindSpeeds_ResultantVector = np.round(np_HourlyWindSpeeds_ResultantVector, decimals=3)
                    np_HourlyWindSpeeds_UnCorrected_ResultantVector = np.round(np_HourlyWindSpeeds_UnCorrected_ResultantVector, decimals=3)

                    np_allMSR_HourlyCF[MSR_Counter] = np_HourlyCF
                    np_allMSR_HourlyWindSpeeds_corrected_ResultantVector[MSR_Counter] = np_HourlyWindSpeeds_ResultantVector
                    np_allMSR_HourlyWindSpeeds_uncorrected_ResultantVector[MSR_Counter] = np_HourlyWindSpeeds_UnCorrected_ResultantVector
                    np_allMSR_ERA5_AnnualMean[MSR_Counter] = ERA5_AnnualMean
                    np_allMSR_Wind_GWA_MSR_Mean[MSR_Counter]= Wind_GWA_MSR_Mean
                print ("%s MSR %s"%(MSR_Counter,RE))

            #prepare the data set and export to xlsx
            pd_output = gpd_MSR_Attributes
            pd_output.rename(columns={"FID": "MSR_ID"}, inplace=True)
            HourTags = ['H%s' % i for i in range(1, 8761)]
            if RE==SolarPVNameConvention:
                print("creating solarpv files")
                pd_SolarBiasInformation = pd.DataFrame(
                    {"ERA_GHI KWh/m2/yr": np_allMSR_GHI_ERA5OriginalAnnualYield, "GSA_GHI KWh/m2/yr": np_allMSR_GHI_GSA_MSR_Mean,
                     'BiasCorrection Adder Wh for solar hours': np_allMSR_BiasCorrection_GHI_Adder_Wh})
                if flag_Diagnosis:
                    pd_output_diagnosis = pd.concat(
                        [pd_output, pd_SolarBiasInformation, pd.DataFrame(np_allMSR_HourlyGHI_corrected_Wh, columns=HourTags)], axis=1)
                    pd_output_diagnosis.to_csv(
                        OutputCountryFolder_UTC+r"\%s %s BiasCorrected ResourceProfiles.csv" % (country,RE),index=False)
                    print(r"\%s %s BiasCorrected ResourceProfiles.csv created" % (country,RE))

                pd_output = pd.concat([pd_output, pd_SolarBiasInformation, pd.DataFrame(np_allMSR_HourlyCF, columns=HourTags)], axis=1)
                FileAddressCountryProfile_UTC=OutputCountryFolder_UTC + r"\%s %s CFs.csv" % (country, RE)
                pd_output.to_csv(FileAddressCountryProfile_UTC, index=False)
                print(r"UTC-->%s %s CFs.csv created" %(country, RE))

            if RE==WindNameConvention:
                print("creating wind files")
                pd_WindBiasInformation = pd.DataFrame(
                    {"GWA Annual MSR Mean m/s": np_allMSR_Wind_GWA_MSR_Mean, "ERA-Raw Annual Mean Speed m/s":np_allMSR_ERA5_AnnualMean})
                if flag_Diagnosis:
                    pd_output_diagnosis = pd.concat(
                        [pd_output, pd_WindBiasInformation, pd.DataFrame(np_allMSR_HourlyWindSpeeds_corrected_ResultantVector, columns=HourTags)], axis=1)
                    pd_output_diagnosis.to_csv(
                        OutputCountryFolder_UTC+r"\%s %s %sm BiasCorrected ResourceProfiles.csv" % (country, RE, WindTurbineHeight_meters), index=False)
                    print(r"\%s %s corrected resource profiles.csv created" % (country,RE))

                    pd_output_diagnosis = pd.concat(
                        [pd_output, pd_WindBiasInformation, pd.DataFrame(np_allMSR_HourlyWindSpeeds_uncorrected_ResultantVector, columns=HourTags)], axis=1)
                    pd_output_diagnosis.to_csv(
                        OutputCountryFolder_UTC + r"\%s %s %sm UnCorrected ResourceProfiles.csv" % (
                        country, RE, WindTurbineHeight_meters), index=False)
                    print(r"\%s %s uncorrected resource profiles.csv created" % (country, RE))

                pd_output = pd.concat([pd_output, pd_WindBiasInformation, pd.DataFrame(np_allMSR_HourlyCF, columns=HourTags)], axis=1)
                FileAddressCountryProfile_UTC=OutputCountryFolder_UTC+r"\%s %s %sm CFs.csv"%(country, RE, WindTurbineHeight_meters)
                pd_output.to_csv(FileAddressCountryProfile_UTC, index=False)
                print(r"UTC-->%s %s CFs.csv created" %(country, RE))

            pd_UTC = pd.read_csv(FileAddressCountryProfile_UTC)

            pd_LocalTime = CreateLocalTimeProfile(pd_UTC, pd_CountryUTC_offsets, country_withspaces)

            pd_LocalTime.to_csv(OutputCountryFolder_LocalTime + "\\" + FileAddressCountryProfile_UTC.rsplit("\\", 1)[1])
            print(r"LocalTime-->%s" % FileAddressCountryProfile_UTC.rsplit("\\", 1)[1])

        except Exception as e:
            print (e)
            traceback.print_exc()
            pd_LogFile = pd_LogFile.append(pd.DataFrame(["%s: Skipped %s profiles" % (country, RE)], columns=['Log']))
            print("%s: Skipped %s profiles" % (country, RE))
            pd_LogFile.to_csv(OutputFolder_UTCProfiles + '\\' + DateTimeStamp + 'ProfileGenerator_LogFile.csv')