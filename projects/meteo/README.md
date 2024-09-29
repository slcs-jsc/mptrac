
# Meteorological Data Retrieval Scripts

This directory contains a set of scripts for retrieving meteorological data (ERA5, GFS, MERRA-2, and NCEP) from different data centers and converting them for use with MPTRAC. Each script is designed to retrieve data from a specific center and perform necessary conversions.

## Scripts Overview

### `get_era5.sh`

- **Description**: This script retrieves [ERA5 reanalysis data](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5)  from the [ECMWF's Climate Data Store (CDS)](https://cds.climate.copernicus.eu/). ERA5 provides a detailed record of the atmosphere's state on a global scale.

- **Usage**: 
  ```bash
  ./get_era5.sh <year> <month> <day> <dir> <res> <dt>
  ```
  - `year`, `month`, `day`: Specify the date for data retrieval.
  - `dir`: Directory to store the retrieved data.
  - `res`: Grid resolution in degrees.
  - `dt`: Timestep in hours.
  
- **Requirements**: Requires the [CDS API](https://cds.climate.copernicus.eu/how-to-api) and Climate Data Operators (CDO) for data access and conversion.

### `get_gfs.sh`

- **Description**: This script retrieves [GFS (Global Forecast System) data](https://www.emc.ncep.noaa.gov/emc/pages/numerical_forecast_systems/gfs.php), a global weather prediction model, and converts it for MPTRAC use.

- **Usage**: 
  ```bash
  ./get_gfs.sh <year> <dir>
  ```
  - `year`: Specify the year for data retrieval.
  - `dir`: Directory to store the processed data.
  
- **Functionality**: Downloads multiple weather parameters and integrates them over different pressure levels for analysis.

### `get_merra2.sh`

- **Description**: Retrieves [MERRA-2 data](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/) from NASA's GES DISC (Goddard Earth Sciences Data and Information Services Center), which provides historical atmospheric reanalysis data.

- **Usage**: 
  ```bash
  ./get_merra2.sh <year> <dir>
  ```
  - `year`: Year for which the data is being retrieved.
  - `dir`: Directory to store the retrieved files.
  
- **Functionality**: This script handles 3D and surface data, processes time steps, and formats it for use in MPTRAC simulations.

### `get_ncep.sh`

- **Description**: Retrieves [NCEP-NCAR (National Centers for Environmental Prediction) Reanalysis 1 data](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html) and converts it for MPTRAC, using reanalysis data on different atmospheric levels.

- **Usage**: 
  ```bash
  ./get_ncep.sh <year> <dir>
  ```
  - `year`: Specify the year for data retrieval.
  - `dir`: Directory to store the processed data.
  
- **Functionality**: Processes different atmospheric parameters, including air temperature, wind speeds, surface pressure, and geopotential height. The script performs vertical interpolation and data merging for MPTRAC compatibility.

## Dependencies

- **CDO (Climate Data Operators)**: Used for converting NetCDF data files.
- **NCO (NetCDF Operators)**: For handling NetCDF file operations.
- **Python**: Required for interacting with data APIs, such as the CDS API for ERA5.
- **CDS API**: For retrieving data from ECMWF's Climate Data Store.

## Notes

- Ensure that you have the necessary access and API keys set up for downloading data from each respective data center.
- The scripts make use of temporary directories for intermediate processing, which are cleaned up automatically after execution.
