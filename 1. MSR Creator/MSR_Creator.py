#BE SURE TO INSTALL THESE LIBRARIES ON THE SERVER / COMPUTER BEFORE ATTEMPTING TO RUN THE CODE

# import relevant libraries

import math
import os
import json
import scipy.stats
import struct, time
from shapely.geometry import LineString, MultiPolygon
from shapely.ops import split
import pandas as pd
import geopandas as gpd
import shutil
import string
import numpy as np
from scipy.ndimage.measurements import label
from pathlib import Path
import rasterio
from rasterio.features import shapes
import xarray
import xrspatial
import richdem as rd
from geocube.api.core import make_geocube
from shapely.geometry import box, mapping
from rasterstats import zonal_stats
from colorama import Fore
from math import ceil


# import warnings
# warnings.filterwarnings("ignore")

def quadrat_cut_geometry(geometry, quadrat_width):
    #Code adopted from OSMNX with alterations:https://osmnx.readthedocs.io/en/stable/index.html

    # create n evenly spaced points between the min and max x and y bounds
    west, south, east, north = geometry.bounds
    x_num = math.floor((east - west) / quadrat_width) + 1
    y_num = math.floor((north - south) / quadrat_width) + 1
    x_points = np.linspace(west, east, num=2+x_num)
    y_points = np.linspace(south, north, num=2+y_num)


    # create a quadrat grid of lines at each of the evenly spaced points
    vertical_lines = [LineString([(x, y_points[0]), (x, y_points[-1])]) for x in x_points]
    horizont_lines = [LineString([(x_points[0], y), (x_points[-1], y)]) for y in y_points]
    lines = vertical_lines + horizont_lines

    # recursively split the geometry by each quadrat line
    for line in lines:
        geometry = MultiPolygon(split(geometry, line))

    return geometry

def PolygonizeResourcePotential(Path_ResourcePotentialRaster,SubFolder_Polygonization, RE_Technology,
               ResourceThreshold, BandCountForMultiResolveAlgorithm,
               MaxAreaToCapMSRs_kM2, MinContiguousAreaSuitableForMSR_kM2):
    # Open the Competitive resource raster for polygonization
    ResourcePotentialRaster = xarray.open_dataarray(Path_ResourcePotentialRaster)
    ResourcePotentialRaster = ResourcePotentialRaster.squeeze("band")

    # Get MaxResourcePixelValue
    np_ResourcePotentialRaster = ResourcePotentialRaster.data * 1
    np_ResourcePotentialRaster[np.isnan(np_ResourcePotentialRaster)] = 0
    MaxResourcePixelValue = np_ResourcePotentialRaster.max()

    Flag_FirstMSR = 1
    for ResourceBand in range(1, BandCountForMultiResolveAlgorithm + 1):
        Path_ResolvedRaster = SubFolder_Polygonization + RE_Technology + 'ResourceBand%s_resolve.tif' % ResourceBand
        Path_SingleBand_InitialMSRs = SubFolder_Polygonization + RE_Technology + 'ResourceBand%s_InitialMSRs.shp' % ResourceBand
        Path_SingleBand_FinalMSRs = SubFolder_Polygonization + RE_Technology + 'ResourceBand%s_FinalMSRs.shp' % ResourceBand

        ResourceBandUpperLimit = ResourceThreshold + (ResourceBand) * (
                    (MaxResourcePixelValue - ResourceThreshold) / BandCountForMultiResolveAlgorithm)
        ResourceBandLowerLimit = ResourceThreshold + (ResourceBand - 1) * (
                    (MaxResourcePixelValue - ResourceThreshold) / BandCountForMultiResolveAlgorithm)

        print("Creating MSRs for band %s : %s to %s " % (
        ResourceBand, ResourceBandLowerLimit, ResourceBandUpperLimit))

        SubsetResourcePotentialRaster = ResourcePotentialRaster * 1
        SubsetResourcePotentialRaster = SubsetResourcePotentialRaster.where(
            ~(ResourcePotentialRaster < ResourceBandLowerLimit), 0)
        SubsetResourcePotentialRaster = SubsetResourcePotentialRaster.where(
            ~(ResourcePotentialRaster > ResourceBandUpperLimit), 0)
        ResolvedRaster = SubsetResourcePotentialRaster.where(~(SubsetResourcePotentialRaster > 0), 1)
        ResolvedRaster.rio.to_raster(Path_ResolvedRaster)

        # below command converts contigious raster pixels to polygon geometries. Transform is crucial input to provide, otherwise the geometries will assume single pixel size as 1 meter
        InitialPolygons = shapes(ResolvedRaster.data.astype('float32'), mask=ResolvedRaster.data == 1,
                                 transform=ResolvedRaster.rio.transform())
        InitialPolygons = ({'properties': {'raster_val': PixelValue}, 'geometry': PolygonGeometry} for
                           counter, (PolygonGeometry, PixelValue) in enumerate(
            InitialPolygons))
        InitialPolygons = list(InitialPolygons)
        time.sleep(1)
        InitialPolygons = gpd.GeoDataFrame.from_features(InitialPolygons, crs="ESRI:54009")

        if not InitialPolygons.empty:
            InitialPolygons = InitialPolygons.explode(ignore_index=True).droplevel(1).reset_index(drop=True)
            InitialPolygons = InitialPolygons.drop(columns=['raster_val'])
            InitialPolygons.to_file(Path_SingleBand_InitialMSRs)

            InitialPolygonsAboveMinAreaThreshold = InitialPolygons[InitialPolygons.area >= MinContiguousAreaSuitableForMSR_kM2 * 1000000]
            if not InitialPolygonsAboveMinAreaThreshold.empty:

                gpd_SingleBand_Final = InitialPolygonsAboveMinAreaThreshold[
                    InitialPolygonsAboveMinAreaThreshold.area <= MaxAreaToCapMSRs_kM2 * 1000000].reset_index(drop=True)
                gpd_SingleBand_ToBeCapped = InitialPolygonsAboveMinAreaThreshold[
                    InitialPolygonsAboveMinAreaThreshold.area > MaxAreaToCapMSRs_kM2 * 1000000].reset_index(drop=True)

                if len(gpd_SingleBand_ToBeCapped) > 0:
                    for i in range(0, len(gpd_SingleBand_ToBeCapped)):
                        print("dividing polygon %s/%s" % (i + 1, len(gpd_SingleBand_ToBeCapped)))
                        gpd_SinglePolygonParts = gpd.GeoDataFrame(crs=gpd_SingleBand_ToBeCapped.crs,
                                                                  geometry=list(
                                                                      quadrat_cut_geometry(
                                                                          gpd_SingleBand_ToBeCapped.geometry.loc[
                                                                              i],
                                                                          np.sqrt(MaxAreaToCapMSRs_kM2) * 1000)))
                        if i == 0 and gpd_SingleBand_Final.empty:
                            gpd_SingleBand_Final = gpd_SinglePolygonParts
                        else:
                            gpd_SingleBand_Final = gpd.overlay(gpd_SinglePolygonParts, gpd_SingleBand_Final,
                                                               how='union')
                            gpd_SingleBand_Final=gpd_SingleBand_Final[gpd_SingleBand_Final.area >= MinContiguousAreaSuitableForMSR_kM2 * 1000000] # To get rid of any small geometries introduced by the above union operation
                gpd_SingleBand_Final[gpd_SingleBand_Final.area >= MinContiguousAreaSuitableForMSR_kM2 * 1000000].to_file(Path_SingleBand_FinalMSRs)

                gpd_SingleBand = gpd_SingleBand_Final[gpd_SingleBand_Final.area >= MinContiguousAreaSuitableForMSR_kM2 * 1000000].reset_index(
                    drop=True)
                if Flag_FirstMSR == 1:
                    gpd_MultiBand = gpd_SingleBand
                    Flag_FirstMSR = 0
                else:
                    gpd_MultiBand = gpd.overlay(gpd_SingleBand, gpd_MultiBand, how='union')
                    gpd_MultiBand = gpd_MultiBand[gpd_MultiBand.area >= MinContiguousAreaSuitableForMSR_kM2 * 1000000]
    try:
        gpd_MultiBand['FID'] = gpd_MultiBand.index
        return gpd_MultiBand
    except:
        print(Fore.RED+"******MSRs not developed*******")
        return 0


def MinimumDistanceOfMSRCentroidFromGivenGeometrySet(gpd_MSRCentroid, gpd_GeometrySet): #Geometry can be point, line or polygon or mixed
    return gpd_GeometrySet.distance(gpd_MSRCentroid).min()

def ComputeLoadCenterAttributesForMSRCentroid(gpd_MSRCentroid, gpd_LoadCenters):
    pd_LoadCenterDistances= gpd_LoadCenters.geometry.distance(gpd_MSRCentroid)/1000
    bestDistance = pd_LoadCenterDistances.min()

    ClosestCityName = gpd_LoadCenters.name_conve[pd_LoadCenterDistances==bestDistance].iloc[0]
    ClosestCityPopulationCount = gpd_LoadCenters.max_pop_al[pd_LoadCenterDistances==bestDistance].iloc[0]
    PopWithin100km = gpd_LoadCenters.max_pop_al[pd_LoadCenterDistances < 100].sum()
    CityCountWithin100km = len(pd_LoadCenterDistances[pd_LoadCenterDistances < 100])

    if CityCountWithin100km <= 10:
        Cities100kM = np.array2string(gpd_LoadCenters.name_conve[pd_LoadCenterDistances < 100].values, separator=",").strip(
            '[]').replace("'", "").replace(",", ", ")
    else:
        Cities100kM = "Above 10 cities"

    return ClosestCityName, ClosestCityPopulationCount,PopWithin100km, CityCountWithin100km,Cities100kM

def RunResourceSufficiencyStage(Path_SuitableAreaResourceRaster,UserResourceThreshold,ResourceLowerLimit,RE_Technology,LandDiscount,
                                cspFootPrintMWpkM2,windRotorDiameter_meters,windSingleTurbineCapacity_Watts,
                                DefaultMinContiguousAreaSuitableForMSR_kM2,CountryArea_kM2):


    #Resource sufficiency check script

    #Reduce cutoff threshold for resource lagging countries until we reach 0 or at yield that meets sufficiency criteria
    SuitableAreaResourceRaster=xarray.open_dataarray(Path_SuitableAreaResourceRaster)
    SuitableAreaResourceRaster=SuitableAreaResourceRaster.squeeze("band")

    RasterPixelSize_m=abs(SuitableAreaResourceRaster.affine[0]) #Raster should be in meter measured Coordiante Reference System (CRS)
    MinContiguousPixelsToRetain = ceil(DefaultMinContiguousAreaSuitableForMSR_kM2*1000000/(RasterPixelSize_m*RasterPixelSize_m))

    np_SuitableAreaResourceRaster = SuitableAreaResourceRaster.to_numpy()
    ResourceThreshold = UserResourceThreshold
    CutoffNormalized=1

    WindLandClasses = np.arange(0, 12.5,0.5)  # Land clusters in m/s [ 0.5,  1. ,  1.5,  2. ,  2.5,  3. ...12.], count=25,actual bins are 24 [0-0.5),[0.5-1)...[11.5-12]
    WindProductionPercentage_perLandClass = np.array([0, 0, 0, 0, 0, 1, 3, 6, 9, 13, 18, 24, 30, 37,43, 48, 54, 58, 61, 64, 65, 66, 66, 65],
                                                     dtype=int)  # % of Nameplate production (Nameplate capacity KW*8760) per 24 land clusters/bins in range[0,12 m/s] at 0.5 step as per IRENA-KTH study:https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2014/IRENA_Africa_Resource_Potential_Aug2014.pdf

    cspLandClasses = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],dtype=int)  # Land clusters in kWh/m2-day , count=18,actual bins are 17
    cspProductionPercentage_perLandClass = np.array([0.00, 7.84, 21.40, 31.01, 38.47, 44.57, 49.72, 54.19, 58.12, 61.65, 64.83, 67.74, 70.42, 72.90,
                                            75.20, 77.36, 79.39, 81.30, ],dtype=float)  # % of Nameplate production (Nameplate capacity KW*8760) per 17 land bins as per IRENA-LBNL study:https://www.irena.org/publications/2015/Oct/Renewable-Energy-Zones-for-the-Africa-Clean-Energy-Corridor

    BreakWhileLoop =0
    #this while loop repeats until sufficiency or lower resource limit reaches.
    while (BreakWhileLoop == 0 and ResourceThreshold>=ResourceLowerLimit):
        # subset resource raster as per resource cutofff
        np_SuitableAreaResourceRaster_Filtered = np.where(np_SuitableAreaResourceRaster < ResourceThreshold, 0,
                                                          np_SuitableAreaResourceRaster)
        # get same shape array with each pixel value replaced by feature no it belongs. Label by default detects features of cross pattern i.e. diagnoal pixels are ignored
        feat, count = label(np_SuitableAreaResourceRaster_Filtered)
        # count no of pixels included per feature. This array is different dimension
        FeaturePixelCount = np.bincount(feat[feat >= 0])
        # exclude features below the minimum area threshold
        DesiredFeatures_toRetain = np.where(FeaturePixelCount > MinContiguousPixelsToRetain)
        # initialize new array of same shape as the filtered resource raster, fill it with original pixel values feature by feature, considering only the retained features.
        np_SuitableAreaResourceRaster_withoutSmallContagiousRegions = np.zeros_like(
            np_SuitableAreaResourceRaster_Filtered)
        if len(DesiredFeatures_toRetain[0]) > 1:
            for f in DesiredFeatures_toRetain[0][1:]:
                np_SuitableAreaResourceRaster_withoutSmallContagiousRegions = np.where(feat == f,np_SuitableAreaResourceRaster_Filtered,np_SuitableAreaResourceRaster_withoutSmallContagiousRegions)
        print("Retained features in total= %s" % len(DesiredFeatures_toRetain[0]))
        np_SuitableAreaResourceRaster_withoutSmallContagiousRegions[np.isnan(np_SuitableAreaResourceRaster_withoutSmallContagiousRegions)] = 0

        # get solarpv, Wind, solar csp yields and compare with user provided requirement
        if RE_Technology == 'solarpv':
            IndicativeYield_GWh = np_SuitableAreaResourceRaster_withoutSmallContagiousRegions.sum() *RasterPixelSize_m*RasterPixelSize_m* (365 / 1000000) * (0.16 / 5) * LandDiscount  #  10e6 for convertion to GWh annual yield, using 16% solar PV conversion efficiency and 20% spacing factor as per KTH report

        if RE_Technology == 'solarcsp':
            #area in kM2
            Area_perSpatialCluster=np.histogram(np_SuitableAreaResourceRaster_withoutSmallContagiousRegions, bins=cspLandClasses)[0] *RasterPixelSize_m*RasterPixelSize_m/1000000
            #Deployable max MWs per spatial cluster assuming LandDiscount
            cspMaxCapacity_perSpatialCluster = cspFootPrintMWpkM2* Area_perSpatialCluster*LandDiscount
            #Get CSP yield in GWh
            IndicativeYield_GWh=(cspMaxCapacity_perSpatialCluster*8760*cspProductionPercentage_perLandClass/100).sum()/1000

        if RE_Technology == 'wind':
            #area in kM2
            Area_perSpatialCluster=np.histogram(np_SuitableAreaResourceRaster_withoutSmallContagiousRegions, bins=WindLandClasses)[0]*RasterPixelSize_m*RasterPixelSize_m/1000000
            #Deployable max MWs per spatial cluster assuming LandDiscount
            WindMaxCapacity_perSpatialCluster = np.round(Area_perSpatialCluster / (5 * windRotorDiameter_meters * 3 * windRotorDiameter_meters / 1000000),0)*(windSingleTurbineCapacity_Watts / 1e6)*LandDiscount
            #Get Wind yield in GWh
            IndicativeYield_GWh=(WindMaxCapacity_perSpatialCluster*8760*WindProductionPercentage_perLandClass/100).sum()/1000
            #print("Area_perSpatialCluster= %s"%Area_perSpatialCluster)
            #print("Wind Land Bins= %s" % WindLandClasses)
            #print("WindMaxCapacity_perSpatialCluster= %s" % WindMaxCapacity_perSpatialCluster)
            #print("production efficiency per Spatial Cluster=%s"%(WindProductionPercentage_perLandClass/100))
            #print("Max production per Spatial Cluster=%s" % (WindMaxCapacity_perSpatialCluster*8760/1000))
            #print("Actual production per Spatial Cluster=%s" % (WindMaxCapacity_perSpatialCluster * 8760*WindProductionPercentage_perLandClass / 100000))
            # mean annual wind speed across all pixels per bin (optional)
            #SpatialAverage_ofAnnualMeanWindSpeed_perSpatialCluster = scipy.stats.binned_statistic(np_SuitableAreaResourceRaster_resolved,
            #							 np_SuitableAreaResourceRaster_resolved, statistic='mean',bins=WindLandClasses)[0]
            #print("SpatialAverage_ofAnnualMeanWindSpeed_perSpatialCluster= %s" % SpatialAverage_ofAnnualMeanWindSpeed_perSpatialCluster)


        print("Indicative %s Yield is %s GWH at resource threshold of %s KWh/m2-d and cutoff of %s" % (RE_Technology,IndicativeYield_GWh,  ResourceThreshold,CutoffNormalized))

        SufficiencyParameter=np.count_nonzero(np_SuitableAreaResourceRaster_withoutSmallContagiousRegions)*RasterPixelSize_m*RasterPixelSize_m/1000000
        SufficiencyCondition=CountryArea_kM2*0.05

        if SufficiencyParameter > SufficiencyCondition:
            BreakWhileLoop = 1
        else:
            ResourceThreshold = ResourceThreshold - 0.01
            if ResourceThreshold>=ResourceLowerLimit:
                CutoffNormalized = (ResourceThreshold- ResourceLowerLimit)/(UserResourceThreshold - ResourceLowerLimit)

    return ResourceThreshold, IndicativeYield_GWh


'''
*********Main Code**********
For detail on methodology, please consult the paper: Link soon coming 

- Script 1: Methodology Stage 1 (part-i) clipping and distance surfaces, Stage 1 (part-ii) Scoring
- Script 2: Methodology Stage 2 get Competitive resource potential with optional Stage-3 resource sufficiency check
- Script 3: Methodology Stage 4 part(i) Polygonization
- Script 4: Methodology Stage 4 part(ii) Attribution


'''
#Read control input file
ControlDataSetNames=pd.read_excel('ControlFile_MSRCreator.xlsx', sheet_name="input dataset names", index_col=0)
ControlCountryWiseInputs=pd.read_excel('ControlFile_MSRCreator.xlsx', sheet_name="country wise inputs", index_col=0)
ControlConfigurations=pd.read_excel('ControlFile_MSRCreator.xlsx', sheet_name="Configurations", index_col=0)
ControlPaths=pd.read_excel('ControlFile_MSRCreator.xlsx', sheet_name="Paths", index_col=0)
ControlAnalysisInputs=pd.read_excel('ControlFile_MSRCreator.xlsx', sheet_name="AnalysisInputs", index_col=0).transpose().drop(index='Comments')

HomeDirectory=str(ControlPaths.loc["HomeDirectory"][0])
InputSpatialDatasetsFolder = HomeDirectory + ControlPaths.loc["FolderAddress_InputSpatialDatasets"][0]

#Fetch run configuration
AllCountries=pd.read_csv(ControlPaths.loc["FileAddress_CountryNamesList"][0],names=["Ct"])
RE_Technology = ControlConfigurations.loc["RE_Technology"][0]
RoadType = ControlConfigurations.loc["RoadType (Include till this type)"][0]
Flag_RelaxThresholdsForResourceLagingCountries = ControlConfigurations.loc["Flag_RelaxThresholdsForResourceLagingCountries"][0]
Flag_RoadsBufferedSearch = ControlConfigurations.loc["Flag_RoadsBufferedSearch"][0]
Flag_GridBufferedSearch = ControlConfigurations.loc["Flag_GridBufferedSearch"][0]
BandCountForMultiResolveAlgorithm = ControlConfigurations.loc["BandCountForMultiResolveAlgorithm"][0]
DefaultMinCapacitySuitableToCreateMSR_MW= ControlConfigurations.loc["DefaultMinCapacitySuitableToCreateMSR_MW"][0]
RunStatus_ofProcessScripts=[ControlConfigurations.loc["Perform Stage 1 Clipping multi-country datasets, prepare distance surfaces, scoring all data layers"][0],
                          ControlConfigurations.loc["Perform Stage 2 Get resource potential with or without resource sufficiency check (stage-3)"][0],
                          ControlConfigurations.loc["Perform Stage 4 part(i) Polygonization"][0],
                          ControlConfigurations.loc["Perform Stage 4 part(ii) Attribution (MSR capacity, area, distance to grid, road and others)"][0]]
Scripts=[]
for i in range(1, 5):
    if not RunStatus_ofProcessScripts[i-1]==0:
        Scripts.append(i)


#assign file names
FileName_PopulationDensity=ControlDataSetNames.loc["FileName_PopulationDensity"][0]
FileName_LandCover=ControlDataSetNames.loc["FileName_LandCover"][0]
FileName_Elevation=ControlDataSetNames.loc["FileName_Elevation"][0]
FileName_ProtectedAreas = ControlDataSetNames.loc["FileName_ProtectedAreas"][0]
FileName_Substations = ControlDataSetNames.loc["FileName_Substations"][0]
FileName_UrbanAreaLoadCenters=ControlDataSetNames.loc["FileName_UrbanAreaLoadCenters"][0]
FileName_Roads = ControlDataSetNames.loc["FileName_Roads"][0]
FileName_PowerGrid=ControlDataSetNames.loc["FileName_PowerGrid"][0]
FileName_TransmissionGrid=ControlDataSetNames.loc["FileName_TransmissionGrid"][0]
FileName_ContinentDistanceSurface_Tgrid=ControlDataSetNames.loc["FileName_ContinentDistanceSurface_Tgrid"][0]
FileName_DistributionGrid=ControlDataSetNames.loc["FileName_DistributionGrid"][0]
FileName_CountryBoundaries=ControlDataSetNames.loc["FileName_CountryBoundaries"][0]
FileName_GHI_Map=ControlDataSetNames.loc["FileName_GHI_Map"][0]
FileName_DNI_Map = ControlDataSetNames.loc["FileName_DNI_Map"][0]
FileName_WindSpeedMap=ControlDataSetNames.loc["FileName_WindSpeedMap"][0]
FileName_WaterBodies=ControlDataSetNames.loc["FileName_WaterBodies"][0]



#assign RE tech related inputs

PV_SlopeThreshold = [int(ControlAnalysisInputs.PV_SlopeThreshold)]  # maximum percentage used
CSP_SlopeThreshold = [int(ControlAnalysisInputs.CSP_SlopeThreshold)]  # maximum percentage used
Wind_SlopeThreshold = [int(ControlAnalysisInputs.Wind_SlopeThreshold)]  # maximum percentage used
PopulationThreshold = [int(ControlAnalysisInputs.PopulationThreshold)]  # The min and max of population to give respectively 0 and 1
GHI_Thresholds = [float(ControlAnalysisInputs.PV_GHI_LowerLimit), float(ControlAnalysisInputs.PV_GHI_Threshold)]  # The min and max of ghi to be converted to 0 and 1 respectively. In between a linear progression will be used
DNI_Thresholds = [float(ControlAnalysisInputs.CSP_DNI_LowerLimit), float(ControlAnalysisInputs.CSP_DNI_Threshold)]  # The min and max of dni to be converted to 0 and 1 respectively. In between a linear progression will be used
Wind_Threshold = [float(ControlAnalysisInputs.Wind_SpeedLowerLimit),float(ControlAnalysisInputs.Wind_SpeedThreshold)]  # The min and max of wind speed to be converted to 0 and 1 respectively. In between a linear progression will be used
LandDiscount_PV= float(int(ControlAnalysisInputs.SolarPV_LandDiscountFactor)/100)
LandDiscount_CSP= float(int(ControlAnalysisInputs.SolarCSP_LandDiscountFactor)/100)
LandDiscount_Wind = float(int(ControlAnalysisInputs.Wind_LandDiscountFactor)/100)
MaxMSRCapactiy = int(int(ControlAnalysisInputs.MSR_MaxCapacityAllowed))

if RE_Technology == 'solarpv':
    ResourceRasterName=FileName_GHI_Map
    LandDiscount = LandDiscount_PV
    SlopeThreshold = PV_SlopeThreshold
    RE_SpatialFootPrint_MWperkM2=int(ControlAnalysisInputs.PV_FootPrint_MWperkM2)
    MaxAreaToCapMSRs_kM2 = MaxMSRCapactiy / LandDiscount / RE_SpatialFootPrint_MWperkM2
    DefaultMinContiguousAreaSuitableForMSR_kM2= DefaultMinCapacitySuitableToCreateMSR_MW / LandDiscount / RE_SpatialFootPrint_MWperkM2
    ResourceLowerLimit = GHI_Thresholds[0]
    UserResourceThreshold = GHI_Thresholds[1]
    RunInfoColumnHeaders=['Resource Threshold kWh/m2day', 'Yield GWh']
if RE_Technology == 'solarcsp':
    ResourceRasterName=FileName_DNI_Map
    LandDiscount = LandDiscount_CSP
    SlopeThreshold = CSP_SlopeThreshold
    RE_SpatialFootPrint_MWperkM2 = int(ControlAnalysisInputs.CSP_FootPrint_MWperkM2)
    MaxAreaToCapMSRs_kM2 = MaxMSRCapactiy / LandDiscount / RE_SpatialFootPrint_MWperkM2
    DefaultMinContiguousAreaSuitableForMSR_kM2 = DefaultMinCapacitySuitableToCreateMSR_MW / LandDiscount / RE_SpatialFootPrint_MWperkM2
    ResourceLowerLimit = DNI_Thresholds[0]
    UserResourceThreshold = DNI_Thresholds[1]
    RunInfoColumnHeaders=['Resource Threshold kWh/m2day', 'Yield GWh']
if RE_Technology == 'wind':
    ResourceRasterName=FileName_WindSpeedMap
    LandDiscount = LandDiscount_Wind
    SlopeThreshold = Wind_SlopeThreshold
    number_of_turbines_per_kM2 = math.floor((1) / (5 * (int(ControlAnalysisInputs.WindTurbineRotorDiameter_meters) / 1000) * 3 * (int(ControlAnalysisInputs.WindTurbineRotorDiameter_meters) / 1000)))
    Wind_FootPrint_MWperkM2 = number_of_turbines_per_kM2 * int(ControlAnalysisInputs.WindTurbineCapacity_Watts) / 1e6
    RE_SpatialFootPrint_MWperkM2=Wind_FootPrint_MWperkM2
    MaxAreaToCapMSRs_kM2 = MaxMSRCapactiy / LandDiscount / Wind_FootPrint_MWperkM2
    DefaultMinContiguousAreaSuitableForMSR_kM2 = DefaultMinCapacitySuitableToCreateMSR_MW / LandDiscount / RE_SpatialFootPrint_MWperkM2
    ResourceLowerLimit = Wind_Threshold[0]
    UserResourceThreshold = Wind_Threshold[1]
    RunInfoColumnHeaders=['Resource Threshold m/s', 'Yield GWh']

rotor_diammeter, turbine_nameplate_capacity = int(ControlAnalysisInputs.WindTurbineRotorDiameter_meters), int(ControlAnalysisInputs.WindTurbineCapacity_Watts)

gdf_CountryBoundaries=gpd.read_file(InputSpatialDatasetsFolder+FileName_CountryBoundaries+".shp")
SubfolderCountryMapsForClipping=HomeDirectory+r"\RegionBoundaryMaps"

pd_LogFile=pd.DataFrame()
for CountryCounter in range(0,len(AllCountries)):#country wise loop
    RegionName_withSpaces=AllCountries.Ct[CountryCounter]
    RegionName_withoutSpaces = AllCountries.Ct[CountryCounter].replace(" ", "")
    print(Fore.GREEN+"Running MSR script for %s"%RegionName_withoutSpaces)

    if Flag_RoadsBufferedSearch:
        RoadsBufferDistance_meters = int(ControlCountryWiseInputs.loc[RegionName_withoutSpaces][0]) * 1000

    if Flag_GridBufferedSearch:
        GridBufferDistance_meters = ControlCountryWiseInputs.loc[RegionName_withoutSpaces][1] * 1000



    #assign and create folders
    OutputFolder = HomeDirectory +ControlPaths.loc["FolderAddress_OutputFolder"][0] + RegionName_withoutSpaces + "/"
    SubfolderStage1_Clipping = OutputFolder + "Stage1_InputDatasets/"
    if not os.path.isdir(SubfolderStage1_Clipping):
        os.makedirs(SubfolderStage1_Clipping)
    SubfolderStage1_Scoring = OutputFolder + "Stage1_ScoredDatasets/"
    if not os.path.isdir(SubfolderStage1_Scoring):
        os.makedirs(SubfolderStage1_Scoring)
    SubfolderStage2_CompetitiveResource = OutputFolder + "Stage2_CompetitiveResourceArea/"
    if not os.path.isdir(SubfolderStage2_CompetitiveResource):
        os.makedirs(SubfolderStage2_CompetitiveResource)
    SubFolder_Polygonization = OutputFolder + "Stage4_Polygonization/"
    if not os.path.isdir(SubFolder_Polygonization):
        os.makedirs(SubFolder_Polygonization)
    SubFolder_Final_MSRs = OutputFolder + "Stage4_MSR/"
    if not os.path.isdir(SubFolder_Final_MSRs):
        os.makedirs(SubFolder_Final_MSRs)

    #Assign File paths
    Path_FinalMSRs = SubFolder_Final_MSRs + RE_Technology + '_FinalMSRs.shp'


    #prepare country boundary shapefile
    gdf_SingleCountry=gdf_CountryBoundaries[gdf_CountryBoundaries.name == RegionName_withSpaces]
    if not os.path.isdir(SubfolderCountryMapsForClipping):
        os.makedirs(SubfolderCountryMapsForClipping)
    gdf_SingleCountry.to_crs('EPSG:4326').to_file(SubfolderCountryMapsForClipping+'/' + RegionName_withoutSpaces + '.shp')
    CountryArea_kM2 = gdf_SingleCountry.to_crs("ESRI:54009").area.iloc[0] / 1000000

    for Script in Scripts:
        if Script == 1:
            UpperLeftX,LowerRightY, LowerRightX, UpperLeftY=gdf_SingleCountry.total_bounds
            MinX, MinY, MaxX, MaxY = gdf_SingleCountry.total_bounds
            ClipGeometry=json.dumps(mapping(box(UpperLeftX, UpperLeftY, LowerRightX, LowerRightY)))

            print(Fore.BLUE+"Starting Stage 1 (part-i) clipping and distance surfaces")

            for Raster in [FileName_PopulationDensity, FileName_LandCover, FileName_Elevation, ResourceRasterName]:
                InputRasterDataset = xarray.open_dataarray("%s%s.tif" % (InputSpatialDatasetsFolder,Raster))
                ClippedRaster = InputRasterDataset.rio.clip_box(MinX, MinY, MaxX, MaxY)  # clip to envelope, saves time
                del (InputRasterDataset)
                ClippedRaster = ClippedRaster.rio.clip(gdf_SingleCountry.geometry) #clip to country boundaries
                ClippedRaster.rio.to_raster("%s%s_%s_clipped.tif" % (SubfolderStage1_Clipping, RE_Technology, Raster))
                ProjectedRaster=ClippedRaster.rio.reproject("ESRI:54009") #we project all rasters to crs measurable in meters as we will deal with meter based thresholds in contrast to degree based
                ProjectedRaster.rio.to_raster("%s%s_%s_projected.tif" % (SubfolderStage1_Clipping, RE_Technology, Raster))
                print("clipped and projected to ESRI:54009 %s raster dataset" % Raster)
            del(ClippedRaster, ProjectedRaster) #to save memory

            for Vector in [FileName_Roads, FileName_WaterBodies, FileName_PowerGrid, FileName_TransmissionGrid, FileName_DistributionGrid, FileName_ProtectedAreas]:
                gdf_ClippedVector=gpd.read_file("%s%s.shp"%(InputSpatialDatasetsFolder,Vector), bbox=tuple(gdf_SingleCountry.total_bounds)) #bbox (bounding box) used to avoid reading unwanted data. Note that this is does not clip features that extend beyond bbox boundary.

                if not gdf_ClippedVector.empty and Vector==FileName_Roads:
                    gdf_ClippedVector = gdf_ClippedVector[gdf_ClippedVector.GP_RTP <= RoadType]

                if not gdf_ClippedVector.empty:
                    gdf_ClippedVector=gpd.clip(gdf_ClippedVector,gdf_SingleCountry.envelope)
                    gdf_ClippedVector["RasterValue"] = 1
                else: # e.g. there are no grid lines within the region. In this case, get polygon covering all region. This whole polygon will be rasterized to 1.
                    gdf_ClippedVector = gpd.GeoDataFrame({'geometry': gdf_SingleCountry.envelope}, geometry='geometry')
                    gdf_ClippedVector[
                        "RasterValue"] = 1  # this column will be the raster value after rasterization. Here, we assume that in absence of data, all bounding box is allowed for deployment.
                    if Vector in [FileName_ProtectedAreas,
                                  FileName_WaterBodies]:  # These are exclusion datasets so the logic needs to inverted. i.e. in absence of data, exclude nothing.
                        gdf_ClippedVector["RasterValue"] = 0

                gdf_ClippedVector = gpd.clip(gdf_ClippedVector, gdf_SingleCountry.geometry)
                if gdf_ClippedVector.empty: #Check again after clipping to borders that if the dataset is again empty then repeat the procedure of self loading appropriate polygon in order for next scoring step to function properly
                    gdf_ClippedVector = gpd.GeoDataFrame({'geometry': gdf_SingleCountry.geometry},geometry='geometry')
                    gdf_ClippedVector["RasterValue"] = 1
                    if Vector in [FileName_ProtectedAreas,FileName_WaterBodies]:
                        gdf_ClippedVector["RasterValue"] = 0
                gdf_ClippedVector=gdf_ClippedVector.to_crs('EPSG:4326')  # All CRS projections must be ensured to have same crs before rasterizing
                gdf_ClippedVector.to_file("%s%s_%s_clipped.shp"%(SubfolderStage1_Clipping,RE_Technology,Vector))
                print("clipped %s vector dataset" % Vector)
                Rasterized_ClippedVector=make_geocube(gdf_ClippedVector,measurements=["RasterValue"],resolution=(0.0025, -0.0025), geom=ClipGeometry).fillna(0)#highest resolution belongs to resource raster. resolution must be in degrees as the vector has same WGS84 degree based CRS.
                Rasterized_ClippedVector=Rasterized_ClippedVector.rio.clip(gdf_SingleCountry.geometry)
                Rasterized_ClippedVector.rio.reproject("ESRI:54009").RasterValue.rio.to_raster("%s%s_%s_RasterizedProjected.tif"%(SubfolderStage1_Clipping,RE_Technology,Vector))
                print("%s vector dataset is rasterized and projected ESRI:54009" % Vector)
            del (gdf_ClippedVector, Rasterized_ClippedVector)  # to save memory

            ElevationRaster=rd.LoadGDAL("%s%s_%s_projected.tif" % (SubfolderStage1_Clipping,RE_Technology, FileName_Elevation)) #loads rd array. This class is just like numpy with extra info. Till date, no utility to convert it to numpy and then xarray.
            SlopeRaster=rd.TerrainAttribute(ElevationRaster, attrib='slope_percentage')
            rd.SaveGDAL("%s%s_slope_projected.tif" % (SubfolderStage1_Clipping,RE_Technology), SlopeRaster)
            del (ElevationRaster, SlopeRaster)

            for DatasetName in [FileName_PowerGrid, FileName_TransmissionGrid, FileName_DistributionGrid, FileName_Roads]:
                Raster = xarray.open_dataarray("%s%s_%s_RasterizedProjected.tif" % (SubfolderStage1_Clipping,RE_Technology, DatasetName))
                DistanceSurface=xrspatial.proximity(Raster.squeeze('band'), distance_metric="EUCLEADIAN")
                DistanceSurface = DistanceSurface.rio.reproject("EPSG:4326")
                DistanceSurface=DistanceSurface.rio.clip(gdf_SingleCountry.envelope)
                DistanceSurface = DistanceSurface.rio.clip(gdf_SingleCountry.geometry)
                DistanceSurface=DistanceSurface.rio.reproject("ESRI:54009")
                DistanceSurface.rio.to_raster("%s%s_DistanceSurface_%s.tif" % (SubfolderStage1_Clipping,RE_Technology,DatasetName))
                print("%s distance surface created" % DatasetName)

            print("Stage 1 (part-i) clipping and distance surfaces finished")

            print(Fore.BLUE+"Starting Stage 1 (part-ii) Scoring")

            for LayerToScoreName in ["%s_projected"%FileName_PopulationDensity, "%s_projected"%FileName_LandCover, "%s_projected"%FileName_Elevation,
                                 "%s_RasterizedProjected"%FileName_ProtectedAreas, "%s_RasterizedProjected"%FileName_WaterBodies,
                                 "DistanceSurface_%s"%FileName_Roads, "DistanceSurface_%s"%FileName_TransmissionGrid,
                                 "slope_projected",
                                 "%s_projected"%ResourceRasterName]:

                LayerToScore = xarray.open_dataarray("%s%s_%s.tif" % (SubfolderStage1_Clipping, RE_Technology,LayerToScoreName))
                ScoredLayer = LayerToScore * 0

                if LayerToScoreName=="%s_projected"%FileName_LandCover:
                    ScoredLayer = ScoredLayer.where(~LayerToScore.isin([11, 14, 20, 30, 110, 120, 130, 140, 150, 180, 190, 200]), 1)

                if LayerToScoreName=="%s_projected"%FileName_Elevation:
                    ScoredLayer = ScoredLayer.where(~(LayerToScore<2000),1)

                if LayerToScoreName=="%s_projected"%FileName_PopulationDensity:
                    ScoredLayer = ScoredLayer.where(~(LayerToScore<=PopulationThreshold[0]),1)

                if LayerToScoreName=="%s_RasterizedProjected"%FileName_ProtectedAreas:
                    ScoredLayer = ScoredLayer.where(~(LayerToScore==0),1)

                if LayerToScoreName=="%s_RasterizedProjected"%FileName_WaterBodies:
                    ScoredLayer = ScoredLayer.where(~(LayerToScore==0),1)

                if LayerToScoreName=="DistanceSurface_%s"%FileName_Roads:
                    if Flag_RoadsBufferedSearch:
                        ScoredLayer = ScoredLayer.where(~(LayerToScore<=RoadsBufferDistance_meters),1)
                    else:
                        ScoredLayer = ScoredLayer.where(~(LayerToScore>=0),1)  #set all to score 1

                if LayerToScoreName=="DistanceSurface_%s"%FileName_TransmissionGrid:
                    if Flag_GridBufferedSearch:
                        ScoredLayer = ScoredLayer.where(~(LayerToScore <=GridBufferDistance_meters), 1)
                    else:
                        ScoredLayer = ScoredLayer.where(~(LayerToScore >= 0), 1)  # set all to score 1

                if LayerToScoreName == "slope_projected":
                    ScoredLayer = ScoredLayer.where(~(LayerToScore <= SlopeThreshold[0]),1)

                if LayerToScoreName == "%s_projected"%ResourceRasterName:
                    ScoredLayer = ScoredLayer.where(~(LayerToScore < ResourceLowerLimit),-1)
                    ScoredLayer = ScoredLayer.where(~(LayerToScore > UserResourceThreshold), 1)
                    ScoredLayer = ScoredLayer.where(~(ScoredLayer == 0), (LayerToScore-ResourceLowerLimit)/(UserResourceThreshold-ResourceLowerLimit))

                ScoredLayer.rio.to_raster("%s%s_%s_Scored.tif" % (SubfolderStage1_Scoring, RE_Technology,LayerToScoreName))
                print("%s scored" % LayerToScoreName)
            print("Stage 1 (part-ii) Scoring finished")

        #Starting Stage 2 Get Competitive resource (inclusive of optional stage3 resource sufficiency check)
        if Script == 2:
            print(Fore.BLUE+"Starting Stage 2 Get Competitive resource")

            ResourceRaster=xarray.open_dataarray("%s%s_%s_projected.tif" % (SubfolderStage1_Clipping, RE_Technology, ResourceRasterName))

            SuitableAreaRasterWithNoExclusions = xarray.open_dataarray("%s%s_%s_projected_Scored.tif" % (SubfolderStage1_Scoring, RE_Technology, ResourceRasterName))


            SuitableAreaRaster = SuitableAreaRasterWithNoExclusions #linear score 0-1 or -1
            ExclusionCount=0
            for ScoredLayerName in ["%s_projected" % FileName_PopulationDensity, "%s_projected" % FileName_LandCover,
                                 "%s_projected" % FileName_Elevation,
                                 "%s_RasterizedProjected" % FileName_ProtectedAreas, "%s_RasterizedProjected" % FileName_WaterBodies,
                                 "DistanceSurface_%s" % FileName_Roads, "DistanceSurface_%s" % FileName_TransmissionGrid,
                                 "slope_projected"]:
                ScoredLayer= xarray.open_dataarray("%s%s_%s_Scored.tif" % (SubfolderStage1_Scoring, RE_Technology,ScoredLayerName))
                if not (RE_Technology!='wind' and ScoredLayerName=="%s_projected"%FileName_Elevation):
                    SuitableAreaRaster=SuitableAreaRaster*ScoredLayer.reindex({'x':SuitableAreaRaster.x, 'y':SuitableAreaRaster.y}, method="nearest") # two xarrays can multiple only when resolution is same

                    CompetitiveAreaRaster=SuitableAreaRaster*0
                    CompetitiveAreaRaster=CompetitiveAreaRaster.where(~(SuitableAreaRaster>=1),1)

                    SuitableAreaResourceRaster=ResourceRaster.where(~(SuitableAreaRaster<=0),0)
                    CompetitiveAreaResourceRaster = ResourceRaster.where(~(SuitableAreaRaster <1), 0)
                ExclusionCount=ExclusionCount+1
                #CompetitiveAreaResourceRaster.rio.to_raster("%s%s_CompetitiveResource_%sExclusions.tif" % (SubfolderStage2_CompetitiveResource, RE_Technology,ExclusionCount))

            SuitableAreaResourceRaster=SuitableAreaResourceRaster.rio.reproject("EPSG:4326")
            SuitableAreaResourceRaster = SuitableAreaResourceRaster.rio.clip(gdf_SingleCountry.geometry)
            SuitableAreaResourceRaster=SuitableAreaResourceRaster.rio.reproject("ESRI:54009")
            SuitableAreaResourceRaster.rio.to_raster("%s%s_SuitableResource.tif" % (SubfolderStage2_CompetitiveResource, RE_Technology))

            CompetitiveAreaResourceRaster = CompetitiveAreaResourceRaster.rio.reproject("EPSG:4326")
            CompetitiveAreaResourceRaster = CompetitiveAreaResourceRaster.rio.clip(gdf_SingleCountry.geometry)
            CompetitiveAreaResourceRaster = CompetitiveAreaResourceRaster.rio.reproject("ESRI:54009")
            CompetitiveAreaResourceRaster.rio.to_raster("%s%s_CompetitiveResource.tif" % (SubfolderStage2_CompetitiveResource, RE_Technology))

            #Stage 2 (get Competitive resource) finished with resource sufficiency check
            if Flag_RelaxThresholdsForResourceLagingCountries:
                print(Fore.BLUE+"Performing resource sufficiency check (stage-3) before finishing the stage-2")

                Path_SuitableAreaResourceRaster = "%s%s_SuitableResource.tif" % (SubfolderStage2_CompetitiveResource, RE_Technology)
                cspFootPrint_MWpkM2 = float(ControlAnalysisInputs.CSP_FootPrint_MWperkM2)
                windRotorDiameter_meters = float(ControlAnalysisInputs.WindTurbineRotorDiameter_meters)
                windSingleTurbineCapacity_Watts = float(ControlAnalysisInputs.WindTurbineCapacity_Watts)
                ResourceThreshold, IndicativeYield_GWh = RunResourceSufficiencyStage(Path_SuitableAreaResourceRaster, UserResourceThreshold, ResourceLowerLimit,RE_Technology,LandDiscount,
                                                                      cspFootPrint_MWpkM2, windRotorDiameter_meters,windSingleTurbineCapacity_Watts,
                                                                      DefaultMinContiguousAreaSuitableForMSR_kM2,CountryArea_kM2)

                if ResourceThreshold < UserResourceThreshold:
                    CompetitiveAreaResourceRaster = SuitableAreaResourceRaster.where(~(SuitableAreaResourceRaster < ResourceThreshold), 0)
                    CompetitiveAreaResourceRaster = CompetitiveAreaResourceRaster.rio.reproject("EPSG:4326")
                    CompetitiveAreaResourceRaster = CompetitiveAreaResourceRaster.rio.clip(gdf_SingleCountry.geometry)
                    CompetitiveAreaResourceRaster = CompetitiveAreaResourceRaster.rio.reproject("ESRI:54009")
                    CompetitiveAreaResourceRaster.rio.to_raster("%s%s_CompetitiveResource_Relaxed.tif" % (SubfolderStage2_CompetitiveResource, RE_Technology))
                    pd_LogFile=pd_LogFile.append(pd.DataFrame(["%s: Resource threshold reduced to: %s"%(RegionName_withoutSpaces,ResourceThreshold)], columns=['Log']))

                pd.DataFrame({'ResourceThreshold':ResourceThreshold, 'IndicativeYield_GWh':IndicativeYield_GWh,'MinimumMSR_CapacityCriteria': DefaultMinCapacitySuitableToCreateMSR_MW}, index=[0]).to_csv("%s%s_Log_ResourceIdentificationPolygonization.csv" % (SubfolderStage2_CompetitiveResource, RE_Technology))
                print(
                    "Stage 2 (get Competitive resource) finished with resource sufficiency check and ResourceThreshold of %s" % ResourceThreshold)

            #Stage 2 (get Competitive resource) finished without resource sufficiency check
            else:
                ResourceThreshold = UserResourceThreshold
                IndicativeYield_GWh=np.nan
                pd.DataFrame({'ResourceThreshold':ResourceThreshold, 'IndicativeYield_GWh':IndicativeYield_GWh,'MinimumMSR_CapacityCriteria': DefaultMinCapacitySuitableToCreateMSR_MW}, index=[0]).to_csv("%s%s_Log_ResourceIdentificationPolygonization.csv"%(SubfolderStage2_CompetitiveResource,RE_Technology))
                print("Stage 2 (get Competitive resource) finished without resource sufficiency check")
        if Script == 3:

            print(Fore.BLUE+"Starting stage 4 part(i) Polygonization")

            InputsFromStage2=pd.read_csv("%s%s_Log_ResourceIdentificationPolygonization.csv"%(SubfolderStage2_CompetitiveResource,RE_Technology))
            ResourceThreshold = InputsFromStage2.ResourceThreshold.values[0]
            IndicativeYield_GWh=InputsFromStage2.IndicativeYield_GWh.values[0]

            tuple_RunInfo = [(ResourceThreshold,IndicativeYield_GWh)]
            RunInfo = pd.DataFrame(tuple_RunInfo, columns=RunInfoColumnHeaders)
            RunInfo.to_csv(SubFolder_Final_MSRs + RE_Technology + 'Run info.csv', index=False)


            if ResourceThreshold<UserResourceThreshold:
                Path_ResourcePotentialRaster = "%s%s_CompetitiveResource_Relaxed.tif" % (SubfolderStage2_CompetitiveResource, RE_Technology)
            else:
                Path_ResourcePotentialRaster = "%s%s_CompetitiveResource.tif" % (SubfolderStage2_CompetitiveResource, RE_Technology)


            PolygonizationSuccessful=0
            MinContiguousAreaSuitableForMSR_kM2=DefaultMinContiguousAreaSuitableForMSR_kM2
            MinCapacitySuitableToCreateMSR_MW=DefaultMinCapacitySuitableToCreateMSR_MW

            gpd_MSRs= PolygonizeResourcePotential(Path_ResourcePotentialRaster, SubFolder_Polygonization, RE_Technology,
                                 ResourceThreshold, BandCountForMultiResolveAlgorithm,
                                 MaxAreaToCapMSRs_kM2,MinContiguousAreaSuitableForMSR_kM2)

            if type(gpd_MSRs)==int:
                if os.path.isfile(Path_FinalMSRs):
                    os.remove(SubFolder_Final_MSRs + RE_Technology + '_FinalMSRs.shp')
                    os.remove(SubFolder_Final_MSRs + RE_Technology + '_FinalMSRs.shx')
                    os.remove(SubFolder_Final_MSRs + RE_Technology + '_FinalMSRs.prj')
                    os.remove(SubFolder_Final_MSRs + RE_Technology + '_FinalMSRs.cpg')
                    os.remove(SubFolder_Final_MSRs + RE_Technology + '_FinalMSRs.dbf')
                print (Fore.RED+"*******Sufficient resource not found to create any MSRs**************")
                pd_LogFile=pd_LogFile.append(pd.DataFrame(["%s: Sufficient resource not found to create any MSRs" % RegionName_withoutSpaces], columns=['Log']))
                break
            else:
                gpd_MSRs.to_file(Path_FinalMSRs)
                print(Fore.GREEN+"MSR creation done")


        if Script == 4:
            print(Fore.BLUE+"Stage 4 part(ii) Attribution: Calculating attributes from zone centroids from Lines, Roads, Power stations, Loadcenters")

            gpd_MSRs=gpd.read_file(Path_FinalMSRs)
            gpd_MSRs['AreakM2']=gpd_MSRs.geometry.area/1000000
            gpd_MSRs['CapacityMW'] = gpd_MSRs['AreakM2'] * LandDiscount * RE_SpatialFootPrint_MWperkM2

            dc_DistanceToRoadsStats_perMSR = zonal_stats(Path_FinalMSRs,"%s%s_DistanceSurface_%s.tif" % (SubfolderStage1_Clipping, RE_Technology, FileName_Roads),
                                                         stats="count min mean max median sum")
            gpd_MSRs['RoadDist']=pd.DataFrame(dc_DistanceToRoadsStats_perMSR)['mean']/1000
            print("Distances to roads inserted")


            gdf_ClippedVector=gpd.read_file("%s%s.shp"%(InputSpatialDatasetsFolder,FileName_TransmissionGrid), bbox=tuple(gdf_SingleCountry.total_bounds)) #bbox (bounding box) used to avoid reading unwanted data. Note that this is does not clip features that extend beyond bbox boundary.
            if not gdf_ClippedVector.empty:
                gdf_ClippedVector = gpd.clip(gdf_ClippedVector, gdf_SingleCountry.geometry)
            if gdf_ClippedVector.empty: #i.e there is no transmission grid
                dc_DistanceTo_TGrid_Stats_perMSR = zonal_stats(Path_FinalMSRs, "%s%s.tif"%(InputSpatialDatasetsFolder,FileName_ContinentDistanceSurface_Tgrid),
                                                               stats="count min mean max median sum")
                gpd_MSRs['T_Dist_gf'] = pd.DataFrame(dc_DistanceTo_TGrid_Stats_perMSR)['mean'] / 1000
                print("Transmission Grid is absent. Distance from nearest crossborder transmission grid is computed")
            else:
                dc_DistanceTo_TGrid_Stats_perMSR = zonal_stats(Path_FinalMSRs,"%s%s_DistanceSurface_%s.tif" % (SubfolderStage1_Clipping, RE_Technology, FileName_TransmissionGrid),
                                                         stats="count min mean max median sum")

                gpd_MSRs['T_Dist_gf'] = pd.DataFrame(dc_DistanceTo_TGrid_Stats_perMSR)['mean'] / 1000
                print("Distances to Transmission lines inserted")


            dc_DistanceTo_DGrid_Stats_perMSR = zonal_stats(Path_FinalMSRs,"%s%s_DistanceSurface_%s.tif" % (SubfolderStage1_Clipping, RE_Technology, FileName_DistributionGrid),
                                                         stats="count min mean max median sum")
            gpd_MSRs['D_Dist_gf']=pd.DataFrame(dc_DistanceTo_DGrid_Stats_perMSR)['mean']/1000
            print("Distances to Distribution lines inserted")

            gpd_MSRs['TD_Dist_gf'] = gpd_MSRs[['T_Dist_gf', 'D_Dist_gf']].min(axis=1)
            print("Distances to closest grid line (Transmission or Distribution) inserted")

            gpd_Substations = gpd.read_file("%s%s.shp"%(InputSpatialDatasetsFolder, FileName_Substations),bbox=tuple(gdf_SingleCountry.total_bounds)).to_crs("ESRI:54009")
            gpd_MSRs['SubstnDist'] = gpd_MSRs.centroid.apply(MinimumDistanceOfMSRCentroidFromGivenGeometrySet, gpd_GeometrySet=gpd_Substations.centroid)/1000
            print("distance to nearest substation inserted")

            gpd_LoadCenters = gpd.read_file("%s%s.shp"%(InputSpatialDatasetsFolder, FileName_UrbanAreaLoadCenters),bbox=tuple(gdf_SingleCountry.total_bounds)).to_crs("ESRI:54009")
            gpd_MSRs['Load_dst'] = gpd_MSRs.centroid.apply(MinimumDistanceOfMSRCentroidFromGivenGeometrySet, gpd_GeometrySet=gpd_LoadCenters.centroid)/1000
            print("distance to nearest load center inserted")

            #get series of tuples each carrying load center related attributes per MSR
            pdSeries_tupple = gpd_MSRs.centroid.apply(ComputeLoadCenterAttributesForMSRCentroid, gpd_LoadCenters=gpd_LoadCenters)
            pd_LoadCenterRelatedAttributes=pd.DataFrame(pdSeries_tupple.tolist(),columns=['ClosestCityName', 'ClosestCityPopulationCount', 'PopWithin100km',
                                                        'CityCountWithin100km', 'Cities100kM'], index=pdSeries_tupple.index)
            gpd_MSRs['City_name']=pd_LoadCenterRelatedAttributes['ClosestCityName']
            gpd_MSRs['City_Pop'] = pd_LoadCenterRelatedAttributes['ClosestCityPopulationCount']
            gpd_MSRs['CtLst100kM']=pd_LoadCenterRelatedAttributes['Cities100kM']
            gpd_MSRs['CtCnt100kM'] = pd_LoadCenterRelatedAttributes['CityCountWithin100km']
            gpd_MSRs['PopIn100kM']=pd_LoadCenterRelatedAttributes['PopWithin100km']
            print("load center related attributes inserted")
            gpd_MSRs.to_file(Path_FinalMSRs)


            print(Fore.GREEN+"Attribution complete")

            pd_LogFile = pd_LogFile.append(pd.DataFrame(["%s:Attribution completed" %RegionName_withoutSpaces], columns=['Log']))
            DateTimeStamp = time.localtime()
            DateTimeStamp="%s%s%s%s%s%s"%(DateTimeStamp.tm_year,DateTimeStamp.tm_mon,DateTimeStamp.tm_mday,DateTimeStamp.tm_hour,DateTimeStamp.tm_min,DateTimeStamp.tm_sec)
            pd_LogFile.to_csv(HomeDirectory+"//"+DateTimeStamp + RE_Technology + '_LogFile.csv')
    #try:shutil.rmtree(SubfolderStage1_Clipping)
    # #except: pass
    # try:shutil.rmtree(SubfolderStage1_Scoring)
    # except: pass
    # try:shutil.rmtree(SubfolderStage2_CompetitiveResource)
    # except: pass
    # try:shutil.rmtree(SubfolderStage2_CompetitiveResource)
    # except: pass
    # try:shutil.rmtree(SubFolder_Polygonization)
    # except: pass

