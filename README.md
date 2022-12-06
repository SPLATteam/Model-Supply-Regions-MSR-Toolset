<img src="https://user-images.githubusercontent.com/77098670/158803434-7125693b-283c-4d36-8780-2c57ad2015c9.gif" width="400" height="200">

## Code Authors:
**MSR_Creator.py, AttributorCombiner.py, Screener.py**: Bilal Hussain\
**Profile_Generator.py**: Bilal Hussain, Yunshu Li, Sebastian Sterl\
**Clustering.py**: Yunshu Li\
**Methodology contributions:** Asami Miketa, Daniel Russo, Pablo Carvajal, Mohamed Elabbas\
**Contact** : [bhussain@irena.org](mailto:bhussain@irena.org), [ssterl@irena.org](mailto:ssterl@irena.org), [y_li@outlook.com](mailto:y_li@outlook.com)

Feel free to drop quick emails of your interest, queries and motivations to extend the tool individually or in collaboration. We would also like to track your feedback for impact reporting purposes particularly on how the tool can/has contributed value to your work. Please also cite the MSR toolset methodology [paper] (https://doi.org/10.1038/s41597-022-01786-5) in your academic/professional works.     

## Introduction:

This is the public repository of Model Supply Regions (MSR) tool. The tool produces geo-referenced geometries of potential utility scale power supply regions of weather dependent renewable technologies and develops their cost and performance attributes to serve as inputs for energy/power system planning models including generation expansion models (e.g. IRENA [SPLAT Models for Africa](https://irena.org/energytransition/Energy-System-Models-and-Data/System-Planning-Test-Model)). The current version can model Solar PV &amp; Wind-onshore technologies. This tool has been applied on 50 countries in Africa. Results of Africa analysis along with methodological details of the tool and definitions of computed attributes is available in the methodology [paper](https://doi.org/10.1038/s41597-022-01786-5). Final outputs are placed in [Zenodo](https://zenodo.org/record/7014609#.Y0fDXXZBw2w). Inputs for MSR creator are available in this sharepoint [link](https://irena.sharepoint.com/:f:/s/EnergyPlanningCapacityBuilding-ExternalSharing/Em2CRHjokAZImIodaK5JgbIBuhli6wNMomtQ3uCzks-Gxw?e=beuTOo).

## How to use MSR toolset:

MSR toolset comprises of five python scripts that are run sequentially. A high level description of each script and process flow is illustrated below.

<img src="https://user-images.githubusercontent.com/77098670/158806616-53ffbad7-beac-4263-8cb2-0f92a8e0aa46.gif" width="650" height="900">

Each script is enclosed in a separate folder along with its respective control file.xlsx. Profile\_Generator.py requires ERA5 reanalysis data as input, which can be downloaded from the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/#!/home) (CDS, after having created an account and followed the instructions on [installing the CDS API](https://cds.climate.copernicus.eu/api-how-to)) using the ERA5 API script available in the respective folder as a starting point.

Folders are numbered as per order of running scheme. All scripts use a list of input regions (countrynames.csv) available in the repository. Users can edit this list to run the toolset on desired number of input regions. For a single country run, just add single region name only. Take care that the input region names align with the input shapefile on region boundaries and relevant control files.

Before running any script in the toolset, make sure that the .py file is placed with **ControlFile.xlsx** in same directory. Also make sure that **ControlFile.xlsx** holds correct file/folder addresses and code/analysis configurations as per user requirements. 

## Primary Dependencies:

Geopandas = 0.9.0 | geocube = 0.1.0 | Matplotlib = 3.4.3 | numpy = 1.20.0 | netCDF4 = 1.5.7 | pandas = 1.3.2 | pyproj = 3.1.0 | pvlib = 0.8.1 | rasterio = 1.2.6 | richdem = 0.3.4 | rioxarray = 0.9.0 | rasterstats = 0.15.0 | scipy = 1.7.1 | Shapely = 1.7.1 | Tslearn = 0.5.2 | xarray = 0.20.2 | xarray-spatial = 0.3.0 |

Library versions mentioned above serve only for guidance. The scripts should be able to run flexibly on any latest library versions with or without minor updates as necessary.

## Versions:

Version 0.1.0 - March 2022

## License:
See [license](https://github.com/bhussain89/TestRepository/blob/main/LICENSE) here

## Disk space requirements:

**MSR\_Creator.py** requires large disk space to store outputs. It creates region wise folders and sub-folders carrying datasets from different steps executed by the code. A single country folder size can vary significantly from 2GB for small countries like Ivory Coast to 20 GB for large countries like Algeria. For a set of 50 Africa countries, the script produced outputs of total 140GB for one technology. Most outputs produced by the script are relevant for diagnostic purposes only. User can find several commented commands at the script end which can be activated to delete some output folders that are not used by next scripts.

**Profile\_Generator.py** also produces outputs that occupy large diskspace. For 50 country Africa run, it produced around 5GB data for single technology. If user decides to use just one of two time zones, the non-selected time zone data can be deleted which halves the output data size.

All remaining scripts of MSR toolset, produce outputs of size well below 1GB.

## Toolset run time:

Indicative code run times for a single country are as follows:

**MSR\_Creator:** 15 minutes per technology

**Profile\_Generator:** Depends on number of MSRs for which profiles are to be created. On average it is around 1 minute/20 MSR.

**AttributorCombiner, Screener, Cluster:** Less than 5 minutes

Total run time (all scripts) for complete Africa dataset (50 countries) took 36 hours for single technology. This run time can significantly decrease if stricter exclusion and resource criteria is adopted in ControFile\_MSRCreator.xlsx.
