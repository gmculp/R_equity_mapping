# R_equity_mapping
This repository contains functions for generating files for equity maps in Leaflet.  It is a work in progress (in other words, there are likely some bugs) and will be updated reglarly. 

At this point in time, the following functions are available:

* ```download_USCB_TIGER_files.R```: Function to download necessary files from USCB TIGER website using state and county FIPS codes. *Note that if you are using this code from within an organization's firewall, the USCB TIGER website (https://www2.census.gov/geo/tiger) may need to be whitelisted.*

* ```generate_USCB_topojson_file.R```: Function to generate tract, ZCTA or PUMA geographies in topojson format using state and county FIPS codes. Includes optional logical arguments to omit unpopulated contiguous areas and areas that are artifacts of county misassignment (e.g., piers in Brooklyn that are misassigned to Manhattan). Topojson is much lighter than geojson.


Required packages that must be installed to run this code:

* ```data.table```: for handling large data.frames more efficiently

* ```geos```: for performing spatial processes

* ```sf```: for reading in shapefiles 

* ```censusapi```: for reading in population data

* ```geojsonio```: for loading generated topojson file for testing
 
* ```igraph```: for collapsing directional multipart polylines in edges files into single part polylines


          
          
Here is a code sample for generating topojson files for the Bay Area:
```
source("R/download_USCB_TIGER_files.R")
source("R/generate_USCB_topojson_file.R")

###specify place to store USCB TIGER files###
USCB_TIGER.path <- "C:/Leaflet_resources/census_files"

###specify place to store network files###
output.path <- "C:/leaflet_resources/topojson_files"

###specify data table containing state and county FIPS codes###
FIPS_dt <- data.table(state=rep("06",9),county=c("001","055","085","013","075","095","041","081","097"))

###automatically download all necessary files from USCB TIGER website###
###you will only have to do this once for each decennial census year###
###for 2010... default###
download_USCB_TIGER_files(FIPS.dt,USCB_TIGER.path,"2010") 
###for 2020###
download_USCB_TIGER_files(FIPS.dt,USCB_TIGER.path,"2020")

###generate topojson for 2010 census tracts###
generate_USCB_topojson_file(FIPS.dt, USCB_TIGER.path, geo.year="2010", geo.type="tract", omit.unpopulated=TRUE, omit.artifacts=TRUE, output.file_name=file.path(output.path,"USCB_06_tract_2010.topojson"), in_clus=2)

###generate topojson for 2020 PUMA###
generate_USCB_topojson_file(FIPS.dt, USCB_TIGER.path, geo.year="2020", geo.type="puma", omit.unpopulated=TRUE, omit.artifacts=TRUE, output.file_name=file.path(output.path,"USCB_06_PUMA_2020.topojson"), in_clus=2)

