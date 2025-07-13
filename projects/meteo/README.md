# Meteorological Data Retrieval Scripts

This directory contains a set of scripts for retrieving meteorological data (ECMWF, ERA5, GFS, MERRA-2, and NCEP) from different data centers and converting them for use with MPTRAC. Each script is designed to retrieve data from a specific center and perform necessary conversions.

## Quick Start

Ensure all dependencies are installed and scripts are executable:

```bash
sudo apt install wget cdo nco python3-pip
pip install cdsapi
chmod +x get_*.sh
```

Example: Retrieve ERA5 data for July 1, 2023, at 1.0° resolution and 3-hourly intervals:

```bash
./get_era5.sh 2023 07 01 ./data 1.0 3
```

## Scripts Overview

### `get_era5.sh`

- **Description**: Retrieves [ERA5 reanalysis data](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5) from the [ECMWF's Climate Data Store (CDS)](https://cds.climate.copernicus.eu/).
  
- **Usage**:
  ```bash
  ./get_era5.sh <year> <month> <day> <output_dir> <resolution_deg> <timestep_hours>
  ```
  
- **Parameters**:
  - `year`, `month`, `day`: Date to retrieve.
  - `output_dir`: Directory to store the retrieved data.
  - `resolution_deg`: Grid resolution in degrees (e.g., `1.0`).
  - `timestep_hours`: Timestep in hours (e.g., `3`).

- **Functionality**: Retrieves and converts data to NetCDF format for MPTRAC.

### `get_gfs.sh`

- **Description**: Retrieves [GFS (Global Forecast System) data](https://www.emc.ncep.noaa.gov/emc/pages/numerical_forecast_systems/gfs.php) and prepares it for MPTRAC use.

- **Usage**:
  ```bash
  ./get_gfs.sh <year> <month> <day> <output_dir>
  ```
  
- **Parameters**:
  - `year`, `month`, `day`: Date to retrieve.
  - `output_dir`: Directory to store the processed data.
- **Functionality**: Downloads weather parameters, integrates data over pressure levels, and outputs NetCDF files.

### `get_ifs.sh`

- **Description**: Downloads [ECMWF high-resolution forecast data (IFS or AIFS)](https://www.ecmwf.int/en/forecasts/datasets/) from [ECMWF Open Data](https://www.ecmwf.int/en/forecasts/datasets/open-data) hosted on AWS.

- **Usage**:
  ```bash
  ./get_ifs.sh <year> <month> <day> <output_dir> <ifs|aifs-single>
  ```
  
- **Parameters**:
  - `year`, `month`, `day`: Forecast date (use `today` for the latest).
  - `output_dir`: Directory to store NetCDF files.
  - `ifs|aifs-single`: Choose between deterministic forecast (IFS) or AI forecast (AIFS).
  
- **Functionality**:
  - Downloads 00 UTC forecasts at 0.25° resolution.
  - IFS: 3-hourly to 144h, 6-hourly to 360h. AIFS: 6-hourly.
  - Converts GRIB2 to NetCDF using CDO.

### `get_merra2.sh`

- **Description**: Retrieves [MERRA-2 data](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/) from NASA's GES DISC.

- **Usage**:
  ```bash
  ./get_merra2.sh <year> <month> <day> <output_dir>
  ```
  
- **Parameters**:
  - `year`, `month`, `day`: Date to retrieve.
  - `output_dir`: Directory to store the data.
  
- **Functionality**: Retrieves 3D and surface data, processes time steps, and outputs MPTRAC-compatible NetCDF files.

### `get_ncep.sh`

- **Description**: Retrieves [NCEP-NCAR Reanalysis 1 data](https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html).

- **Usage**:
  ```bash
  ./get_ncep.sh <year> <output_dir>
  ```
  
- **Parameters**:
  - `year`: Year to retrieve.
  - `output_dir`: Directory to store the data.
  
- **Functionality**: Processes air temperature, wind, pressure, and geopotential height. Performs vertical interpolation and conversion for MPTRAC.

## Dependencies

- `wget`: For downloading files.
- `CDO` (Climate Data Operators): Converts GRIB/NetCDF formats.
- `NCO` (NetCDF Operators): Handles NetCDF operations.
- `Python`: Required for API interaction.
- `CDS API`: Needed for ERA5 access.

## Notes

- Ensure all required API keys are set up (e.g., `.cdsapirc` for ERA5).
- Scripts create temporary directories during processing, which are removed after execution.
- All output is in NetCDF format.
- Only `get_ifs.sh` supports `today` as an input for the date.
