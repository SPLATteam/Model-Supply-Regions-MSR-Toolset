import pandas as pd
import os
import geopandas as gpd

# import warnings
# warnings.filterwarnings("ignore")

def extract_excel(OutputFileName, SourceFolder_Profiles, gdf_destination, REtechnology, ProfilesFileName, Flag_RunExcelExtraction,SolarMultiple=2.1):
    if Flag_RunExcelExtraction==1:
        countries= gdf_destination.CtryName.unique()
        gdf_destination.rename(columns={"FID":"MSR_ID"} , inplace=True)

        if REtechnology=="Solar PV":
            pd_SplatReady=pd.DataFrame()
            for ctry in countries:
                pd_Profiles_SingleCountry = pd.read_csv("%s\\%s\\%s %s" % (SourceFolder_Profiles, ctry, ctry, ProfilesFileName))
                LongitudeColumn=pd_Profiles_SingleCountry['Longitude']
                LatitudeColumn=pd_Profiles_SingleCountry['Latitude']
                pd_Profiles_SingleCountry=pd_Profiles_SingleCountry.drop(pd_Profiles_SingleCountry.iloc[:,1:20].columns,axis=1)
                pd_Profiles_SingleCountry['Longitude']=LongitudeColumn
                pd_Profiles_SingleCountry['Latitude']=LatitudeColumn
                pd_SplatReady_SingleCountry = pd.merge(gdf_destination[gdf_destination.CtryName==ctry], pd_Profiles_SingleCountry, on="MSR_ID")

                pd_SplatReady_SingleCountry = pd_SplatReady_SingleCountry.drop(["geometry"], axis=1)
                pd_SplatReady= pd_SplatReady.append(pd_SplatReady_SingleCountry)
                print("Appending csv rows for: %s"%ctry)

            cols = pd_SplatReady.columns.tolist()
            cols=cols[:1] + cols[-2:] + cols[1:-2]
            pd_SplatReady=pd_SplatReady[cols]

            pd_SplatReady.to_csv(OutputFileName)

        if REtechnology=="Wind":
            pd_SplatReady=pd.DataFrame()
            for ctry in countries:
                pd_Profiles_SingleCountry = pd.read_csv("%s\\%s\\%s %s" % (SourceFolder_Profiles, ctry, ctry, ProfilesFileName))
                LongitudeColumn=pd_Profiles_SingleCountry['Longitude']
                LatitudeColumn=pd_Profiles_SingleCountry['Latitude']
                pd_Profiles_SingleCountry=pd_Profiles_SingleCountry.drop(pd_Profiles_SingleCountry.iloc[:,1:19].columns,axis=1)
                pd_Profiles_SingleCountry['Longitude']=LongitudeColumn
                pd_Profiles_SingleCountry['Latitude']=LatitudeColumn

                pd_SplatReady_SingleCountry = pd.merge(gdf_destination[gdf_destination.CtryName==ctry], pd_Profiles_SingleCountry, on="MSR_ID")

                pd_SplatReady_SingleCountry = pd_SplatReady_SingleCountry.drop(["geometry"], axis=1)
                pd_SplatReady= pd_SplatReady.append(pd_SplatReady_SingleCountry)
                print("Appending csv rows for: %s"%ctry)

            cols = pd_SplatReady.columns.tolist()
            cols=cols[:1] + cols[-2:] + cols[1:-2]
            pd_SplatReady=pd_SplatReady[cols]

            pd_SplatReady.to_csv(OutputFileName)


ControlPathsAndConfigurations=pd.read_excel('ControlFile_Screener.xlsx', sheet_name="PathsAndConfig", index_col=0)

OutputFolder=ControlPathsAndConfigurations.loc["OutputFolder"][0]
if not os.path.isdir(OutputFolder):
    os.makedirs(OutputFolder)


SolarPVSourceFolderCarryingProfiles = ControlPathsAndConfigurations.loc["SolarPVSourceFolderCarryingProfiles"][0]
WindSourceFolderCarryingProfiles = ControlPathsAndConfigurations.loc["WindSourceFolderCarryingProfiles"][0]
SolarPV_ProfilesFileName=ControlPathsAndConfigurations.loc["SolarPV_ProfilesFileName"][0]
WindCF_ProfilesFileName=ControlPathsAndConfigurations.loc["WindCF_ProfilesFileName"][0]


Countries=pd.read_csv(ControlPathsAndConfigurations.loc["FileAddress_CountryNamesList"][0],names=["Ct"])
RegionBoundariesShapeFile=ControlPathsAndConfigurations.loc["RegionBoundariesShapeFile"][0]


#Pre-Screen file addresses
SolarPV_PreScreenFile=ControlPathsAndConfigurations.loc["SolarPV_PreScreenFile"][0]
WindPreScreenFile=ControlPathsAndConfigurations.loc["WindPreScreenFile"][0]

# All screening is done after sorting MSRs in descending order of LCOE (supply+transmission+road)
# This script is drafted in such a way that it facilitates the further expansion to any no of screening options

#Current screening options
# 1 Country specific covered area cutoff (Select best MSRs that cover x % of country area)

#read screening options and criteria to apply
ScreeningOptions=pd.read_excel('ControlFile_Screener.xlsx', sheet_name="Select screening option", index_col=0)
CountrySpecificCriteriaSolarPV=pd.read_excel('ControlFile_Screener.xlsx', sheet_name="SolarPV country specific", index_col=0)
CountrySpecificCriteriaWind=pd.read_excel('ControlFile_Screener.xlsx', sheet_name="Wind country specific", index_col=0)



#Read which screening options has been chosen by the user to run (1 means run, 0 means dont run)
ScreeningOptionsSelected= ScreeningOptions[ScreeningOptions['Selection Status']==1].index.astype('int').to_list()
print(ScreeningOptionsSelected)
ExcelExtractionOptionsSelected = ScreeningOptions[ScreeningOptions['Extract in Excel'] == 1].index.astype('int').to_list()
print(ExcelExtractionOptionsSelected)

Flag_RunSolarPV=ControlPathsAndConfigurations.loc["Run code for SolarPV"][0]
Flag_RunWind=ControlPathsAndConfigurations.loc["Run code for Wind"][0]
WindHeight=ControlPathsAndConfigurations.loc["wind hub height in meters"][0]
WindAnnualCFColName = 'CF%sm' % WindHeight
WindAnnualYieldColName = 'Y_GWh%sm' % WindHeight

#Solar MSR Screener
if Flag_RunSolarPV:
    REtechnology="Solar PV"
    gdf= gpd.GeoDataFrame()
    file_to_read=SolarPV_PreScreenFile
    gdf_source = gpd.read_file(file_to_read)
    gdf_source.CF=gdf_source.CF.astype(float)

    # sort in descending order of LCOE
    gdf_source=gdf_source.sort_values(by=['LCOE-MWh'], ascending=True, ignore_index=True)

    for method in ScreeningOptionsSelected:
        print("SolarPV Screening Method %s"%method)
        if method in ExcelExtractionOptionsSelected:
            Flag_RunExcelExtraction=1
            print("Extract Excel")
        else:
            Flag_RunExcelExtraction=0
            print("Not extracting Excel")

        if method==1:
            gdf_destination = gpd.GeoDataFrame()

            for i in range(0, len(Countries)):
                country = Countries.Ct[i]
                country_withoutSpace = country.replace(" ", "")

                gdf_RegionBoundaries = gpd.read_file(RegionBoundariesShapeFile)
                CountryArea_kM2=gdf_RegionBoundaries[gdf_RegionBoundaries.name == country].to_crs("ESRI:54009").area.iloc[0] / 1000000
                cutoff=(CountrySpecificCriteriaSolarPV.loc[country][0]/100)*CountryArea_kM2

                gdf_SingleCountry = gdf_source[gdf_source.CtryName == country_withoutSpace]
                gdf_SingleCountry['CumAreakM2'] = gdf_SingleCountry.AreakM2.cumsum()

                gdf_SingleCountry = gdf_SingleCountry[gdf_SingleCountry['CumAreakM2'] <= cutoff]


                gdf_destination=gpd.GeoDataFrame(pd.concat([gdf_destination, gdf_SingleCountry]))
                print (country, CountryArea_kM2)

            gdf_destination=gdf_destination.drop(['CumAreakM2'], axis=1)
            gdf_destination.to_file(OutputFolder + "\\SolarPV_BestMSRsToCover5%CountryArea.shp")
            extract_excel(OutputFolder + "\\SolarPV_BestMSRsToCover5%CountryArea.csv",
            SolarPVSourceFolderCarryingProfiles, gdf_destination, REtechnology, SolarPV_ProfilesFileName,
            Flag_RunExcelExtraction)

#Wind MSR Screener
if Flag_RunWind:
    REtechnology = "Wind"
    gdf= gpd.GeoDataFrame()
    file_to_read=WindPreScreenFile
    gdf_source = gpd.read_file(file_to_read)
    gdf_source[WindAnnualCFColName] = gdf_source[WindAnnualCFColName].astype(float)

    # sort in descending order of LCOE
    gdf_source=gdf_source.sort_values(by=['LCOE-MWh'], ascending=True, ignore_index=True)

    for method in ScreeningOptionsSelected:
        print("Wind Screening Method %s"%method)
        if method in ExcelExtractionOptionsSelected:
            Flag_RunExcelExtraction=1
            print ("Extract Excel")
        else:
            Flag_RunExcelExtraction=0

        if method==1:
            gdf_destination = gpd.GeoDataFrame()

            for i in range(0, len(Countries)):
                country = Countries.Ct[i]
                country_withoutSpace = country.replace(" ", "")

                gdf_RegionBoundaries = gpd.read_file(RegionBoundariesShapeFile)
                CountryArea_kM2=gdf_RegionBoundaries[gdf_RegionBoundaries.name == country].to_crs("ESRI:54009").area.iloc[0] / 1000000
                cutoff=(CountrySpecificCriteriaSolarPV.loc[country][0]/100)*CountryArea_kM2

                gdf_SingleCountry = gdf_source[gdf_source.CtryName == country_withoutSpace]
                gdf_SingleCountry['CumAreakM2'] = gdf_SingleCountry.AreakM2.cumsum()

                gdf_SingleCountry = gdf_SingleCountry[gdf_SingleCountry['CumAreakM2'] <= cutoff]

                gdf_destination=gpd.GeoDataFrame(pd.concat([gdf_destination, gdf_SingleCountry]))
                print (country, CountryArea_kM2)

            gdf_destination=gdf_destination.drop(['CumAreakM2'], axis=1)
            gdf_destination.to_file(OutputFolder + "\\Wind_BestMSRsToCover5%CountryArea.shp")
            extract_excel(OutputFolder + "\\Wind_BestMSRsToCover5%CountryArea.csv",
            WindSourceFolderCarryingProfiles, gdf_destination, REtechnology, WindCF_ProfilesFileName,
            Flag_RunExcelExtraction)