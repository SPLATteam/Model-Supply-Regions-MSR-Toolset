import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from rasterstats import zonal_stats
import os
import math
import glob
import pandas as pd
import geopandas as gpd
import time


#Read control input file
ControlPathsAndNames=pd.read_excel('ControlFile_AttributorCombiner.xlsx', sheet_name="PathsAndNames", index_col=0)
CostParameters=pd.read_excel('ControlFile_AttributorCombiner.xlsx', sheet_name="CostParameters", index_col=0)
ControlConfigurations=pd.read_excel('ControlFile_AttributorCombiner.xlsx', sheet_name="configurations", index_col=0)
Countries=pd.read_csv(ControlPathsAndNames.loc["FileAddress_CountryNamesList"][0],names=["Ct"])
Input_MSR_Folder=ControlPathsAndNames.loc["Input_MSR_Folder"][0]
OutputFolder=ControlPathsAndNames.loc["OutputFolder"][0]
MSR_ShapeFileNameSuffix=ControlPathsAndNames.loc["MSR_ShapeFileNameSuffix"][0]
WindHeight=round(ControlConfigurations.loc["wind hub height in meters"][0])


ProfileGeneratorControlFile_ControlPathsAndNames=pd.read_excel(ControlPathsAndNames.loc["ProfileGeneratorControlFile"][0], sheet_name="PathsAndNames", index_col=0)

SolarPVNameConvention=ProfileGeneratorControlFile_ControlPathsAndNames.loc["SolarPVNameConvention"][0]
WindNameConvention=ProfileGeneratorControlFile_ControlPathsAndNames.loc["WindNameConvention"][0]
ResourceRasterCarryingSubFolderName=ProfileGeneratorControlFile_ControlPathsAndNames.loc["ResourceRasterCarryingSubFolderName"][0]
MSR_DataCarryingSubFolderName=ProfileGeneratorControlFile_ControlPathsAndNames.loc["MSR_DataCarryingSubFolderName"][0]

# set code run configurations
Input_Profiles_Folder = ControlPathsAndNames.loc["Input_Profiles_Folder"][0]
SolarPV_ProfilesFileName = '%s CFs.csv'%SolarPVNameConvention
WindCF_ProfilesFileName = '%s %sm CFs.csv'%(WindNameConvention, WindHeight)

Flag_PrepareMSR_CostCharachteristics= ControlConfigurations.loc["Compute MSR cost attributes"][0]
Flag_RunSolarPV=ControlConfigurations.loc["Run solar pv"][0]
Flag_RunWind=ControlConfigurations.loc["Run wind"][0]


if not os.path.isdir(OutputFolder):
    os.makedirs(OutputFolder)

pd_LogFile=pd.DataFrame()
DateTimeStamp = time.localtime()
DateTimeStamp = "%s%s%s%s%s%s" % (DateTimeStamp.tm_year, DateTimeStamp.tm_mon, DateTimeStamp.tm_mday, DateTimeStamp.tm_hour, DateTimeStamp.tm_min,DateTimeStamp.tm_sec)

if Flag_RunSolarPV:
    gdf= gpd.GeoDataFrame()
    for i in range (0,len(Countries)):
        country_WithSpaces=Countries.Ct[i]
        country = country_WithSpaces.replace(" ", "")
        try:
            if glob.glob("%s\\%s\\%s\\%s"%(Input_MSR_Folder, country, MSR_DataCarryingSubFolderName,SolarPVNameConvention) + MSR_ShapeFileNameSuffix): #single iteration loop
                file_to_read=glob.glob("%s\\%s\\%s\\%s"%(Input_MSR_Folder, country, MSR_DataCarryingSubFolderName,SolarPVNameConvention) + MSR_ShapeFileNameSuffix)[0]
                print("SolarPV:"+file_to_read)
                gdf_SingleCountry = gpd.read_file(file_to_read)
                gdf_SingleCountry=gdf_SingleCountry.sort_values(by=['FID'])

                gdf_SingleCountry['CtryName'] = country

                dc_GHI_ZoneStats_withSolarGIS = zonal_stats(file_to_read,r"%s\%s\%s\%s_GHI_projected.tif"%(Input_MSR_Folder,country,ResourceRasterCarryingSubFolderName,SolarPVNameConvention),stats="count min mean max median sum")

                MeanResource_AverageAcrossMSR=pd.DataFrame(index=gdf_SingleCountry.index,columns=["ZoneMean_KWh/d-m2"])
                for j in range (0, len(gdf_SingleCountry)):
                    MeanResource_AverageAcrossMSR.loc[j] = dc_GHI_ZoneStats_withSolarGIS[j]['mean']
                gdf_SingleCountry['GHIkWhm2d'] = MeanResource_AverageAcrossMSR


                pd_SolarProfiles=pd.read_csv("%s\\%s\\%s %s"%(Input_Profiles_Folder, country,country,SolarPV_ProfilesFileName))
                gdf_SingleCountry['RawERA_GHI'] = pd_SolarProfiles['ERA_GHI KWh/m2/yr']/365 #daily mean
                gdf_SingleCountry['CorAdderWh'] = pd_SolarProfiles['BiasCorrection Adder Wh for solar hours']
                gdf_SingleCountry['CF']=0.92*pd_SolarProfiles.iloc[:,-8760:].sum(axis=1)*100/8760 #ensure excel table is FID sorted, 4% outage and 4% inverter and wiring losses applied as per IRENA-LBNL report
                gdf_SingleCountry['Y_GWh']=gdf_SingleCountry['CapacityMW']*gdf_SingleCountry['CF']*8760/100000


                if Flag_PrepareMSR_CostCharachteristics:
                    #LCOEs calculated per kWH and scaled to MWh
                    #first element is in kWH so must be multiplied with 1000 to get MW basis. 2nd element is MW basis, so simply added.
                    gdf_SingleCountry['sLCOE-MWh'] = 1000 * ((CostParameters.loc["SolarPV"][0] * CostParameters.loc["SolarPV"][3] + CostParameters.loc["SolarPV"][1]/1000) / (8760 * gdf_SingleCountry['CF'] / 100)) + CostParameters.loc["SolarPV"][2]
                    #TL parameters are in MW basis so the LCOE inherently comes in MWh. Substation costs are for 2 substations suitable for an average 50MW POA/plant. So they must be converted to per MW basis.
                    gdf_SingleCountry['tLCOE-MWh'] = ((CostParameters.loc["SolarPV"][5]* CostParameters.loc["SolarPV"][8]+CostParameters.loc["SolarPV"][6]) * gdf_SingleCountry['T_Dist_gf']+CostParameters.loc["SolarPV"][7]*CostParameters.loc["SolarPV"][8]) / (8760 * gdf_SingleCountry['CF'] / 100)
                    gdf_SingleCountry['tCAPEX-kW'] = (CostParameters.loc["SolarPV"][5]*gdf_SingleCountry['T_Dist_gf'] +
                                                      CostParameters.loc["SolarPV"][7]) / 1000
                    #road costs are for 50MW plant i.e. 50*8760*CF MWh
                    gdf_SingleCountry['rLCOE-MWh'] = ((CostParameters.loc["SolarPV"][10]* CostParameters.loc["SolarPV"][12]+CostParameters.loc["SolarPV"][11])*gdf_SingleCountry['RoadDist']) / (8760 *50* gdf_SingleCountry['CF'] / 100)
                    gdf_SingleCountry['rCAPEX-kW'] = (CostParameters.loc["SolarPV"][10] * gdf_SingleCountry['RoadDist']) /50000

                    gdf_SingleCountry['LCOE-MWh']= gdf_SingleCountry['tLCOE-MWh']+ gdf_SingleCountry['rLCOE-MWh']+gdf_SingleCountry['sLCOE-MWh']
                    gdf_SingleCountry['trCAPEX-kW'] = gdf_SingleCountry['rCAPEX-kW']+gdf_SingleCountry['tCAPEX-kW']


                gdf = gpd.GeoDataFrame(pd.concat([gdf, gdf_SingleCountry]))
        except:
            print("Skipped %s" % country)
            pd_LogFile=pd_LogFile.append(pd.DataFrame(["%s: Skipped %s" % (SolarPVNameConvention, country)], columns=['Log']))
            pass
    gdf.to_file(OutputFolder+"\\SolarPV_prescreen.shp")
    pd_LogFile.to_csv(OutputFolder + '\\'+DateTimeStamp+'solarpv_LogFile.csv')

if Flag_RunWind:
    gdf= gpd.GeoDataFrame()
    for i in range (0,len(Countries)):
        country_WithSpaces=Countries.Ct[i]
        country = country_WithSpaces.replace(" ", "")

        try:
            if glob.glob("%s\\%s\\%s\\%s"%(Input_MSR_Folder, country, MSR_DataCarryingSubFolderName, WindNameConvention) + MSR_ShapeFileNameSuffix): #single iteration loop
                file_to_read=glob.glob("%s\\%s\\%s\\%s"%(Input_MSR_Folder, country, MSR_DataCarryingSubFolderName, WindNameConvention) + MSR_ShapeFileNameSuffix)[0]
                print("Wind:"+file_to_read)
                gdf_SingleCountry = gpd.read_file(file_to_read)
                gdf_SingleCountry=gdf_SingleCountry.sort_values(by=['FID'])

                gdf_SingleCountry['CtryName'] = country

                dc_WS_ZoneStats_withGWA = zonal_stats(file_to_read,r"%s\%s\%s\%s_GWA_Africa100m_projected.tif"%(Input_MSR_Folder,country,ResourceRasterCarryingSubFolderName,WindNameConvention),stats="count min mean max median sum")

                MeanResource_AverageAcrossMSR=pd.DataFrame(index=gdf_SingleCountry.index,columns=["ZoneMean_m/s"])
                #MeanResource_AverageAcrossMSR['ZoneMean_m/s'][pd.isnull(MeanResource_AverageAcrossMSR['ZoneMean_m/s'])] = MeanResource_AverageAcrossMSR['ZoneMean_m/s'].mean() # there can be some very rare cases Nonetypes get involved from zonal_stats command

                for j in range (0, len(gdf_SingleCountry)):
                    MeanResource_AverageAcrossMSR.loc[j] = dc_WS_ZoneStats_withGWA[j]['mean']
                gdf_SingleCountry['MeanSpeed'] = MeanResource_AverageAcrossMSR

                gdf_SingleCountry['IEC_Class'] = "" #just for initialization
                gdf_SingleCountry.loc[gdf_SingleCountry['MeanSpeed']<=7.5,'IEC_Class']="Class-3"
                gdf_SingleCountry.loc[gdf_SingleCountry['MeanSpeed'] >= 8.5,'IEC_Class']= "Class-1"
                gdf_SingleCountry.loc[(gdf_SingleCountry['MeanSpeed'] > 7.5)&(gdf_SingleCountry['MeanSpeed'] < 8.5),'IEC_Class'] = "Class-2"

                AnnualCFColName = 'CF%sm' % WindHeight
                AnnualYieldColName = 'Y_GWh%sm' % WindHeight

                pd_WindProfiles=pd.read_csv("%s\\%s\\%s %s"%(Input_Profiles_Folder, country,country,WindCF_ProfilesFileName))
                gdf_SingleCountry['ERA_WSpeed'] = pd_WindProfiles['ERA-Raw Annual Mean Speed m/s']
                gdf_SingleCountry[AnnualCFColName]=0.83*pd_WindProfiles.iloc[:,-8760:].sum(axis=1)*100/8760 #ensure excel table is FID sorted, 2% outage loss and 15% array and collection losses applied as per IRENA-LBNL report
                gdf_SingleCountry[AnnualYieldColName]=gdf_SingleCountry['CapacityMW']*gdf_SingleCountry[AnnualCFColName]*8760/100000


                if Flag_PrepareMSR_CostCharachteristics:
                    var1=["Class-3","Class-2","Class-1"]
                    var2=["Wind_Class3","Wind_Class2","Wind_Class1"]
                    for i in [0,1,2]:
                        #LCOEs calculated per kWH and scaled to MWh
                        #first element is in kWH so must be multiplied with 1000 to get MW basis. 2nd element is MW basis, so simply added.
                        gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class']==var1[i],'sLCOE-MWh'] = 1000 * ((CostParameters.loc[var2[i]][0] * CostParameters.loc[var2[i]][3] + CostParameters.loc[var2[i]][1]/1000) / (8760 * gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class']==var1[i],AnnualCFColName] / 100)) + CostParameters.loc[var2[i]][2]
                        #TL parameters are in MW basis so the LCOE inherently comes in MWh. Substation costs are for 2 substations suitable for an average 50MW POA/plant. So they must be converted to per MW basis.
                        gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class']==var1[i],'tLCOE-MWh'] = ((CostParameters.loc[var2[i]][5]* CostParameters.loc[var2[i]][8]+CostParameters.loc[var2[i]][6]) * gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class']==var1[i],'T_Dist_gf']+CostParameters.loc[var2[i]][7]*CostParameters.loc[var2[i]][8]) / (8760 * gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class']==var1[i],AnnualCFColName] / 100)
                        gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class'] == var1[i], 'tCAPEX-kW'] = (CostParameters.loc[var2[i]][5]
                                                                                                         *gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class'] ==var1[i], 'T_Dist_gf'] +
                                                                                                         CostParameters.loc[var2[i]][7]) / 1000

                        #road costs are for 50MW plant i.e. 50*8760*CF MWh
                        gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class']==var1[i],'rLCOE-MWh'] = ((CostParameters.loc[var2[i]][10]* CostParameters.loc[var2[i]][12]+CostParameters.loc[var2[i]][11])*gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class']==var1[i],'RoadDist']) / (8760 *50* gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class']==var1[i],AnnualCFColName] / 100)
                        gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class']==var1[i],'rCAPEX-kW'] = (CostParameters.loc[var2[i]][10]*gdf_SingleCountry.loc[gdf_SingleCountry['IEC_Class']==var1[i],'RoadDist']) /50000


                    gdf_SingleCountry['LCOE-MWh']= gdf_SingleCountry['tLCOE-MWh']+ gdf_SingleCountry['rLCOE-MWh']+gdf_SingleCountry['sLCOE-MWh']
                    gdf_SingleCountry['trCAPEX-kW'] = gdf_SingleCountry['rCAPEX-kW'] + gdf_SingleCountry['tCAPEX-kW']

                gdf = gpd.GeoDataFrame(pd.concat([gdf, gdf_SingleCountry]))

        except:
            print ("Skipped %s"%country)
            pd_LogFile=pd_LogFile.append(pd.DataFrame(["%s: Skipped %s" % (WindNameConvention,country)],columns=['Log']))
            pass

    gdf.to_file(OutputFolder+"\\Wind_prescreen.shp")

    pd_LogFile.to_csv(OutputFolder + '\\'+DateTimeStamp+'Wind_LogFile.csv')