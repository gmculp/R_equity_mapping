# R_equity_mapping
This repository contains functions for generating files for equity maps in Leaflet.  It is a work in progress (in other words, there are likely some bugs) and will be updated reglarly. 

At this point in time, the following functions are available:

* ```download_USCB_TIGER_files.R```: Function to download necessary files from USCB TIGER website using state and county FIPS codes. *Note that if you are using this code from within an organization's firewall, the USCB TIGER website (https://www2.census.gov/geo/tiger) may need to be whitelisted.*

* ```generate_USCB_topojson_file.R```: Function to generate tract, ZCTA or PUMA geographies in topojson format. Includes optional logical arguments to omit unpopulated contiguous areas and areas that are artifacts of county misassignment (e.g., piers in Brooklyn that are misassigned to Manhattan). Topojson is much lighter than geojson.
