library(terra)
library(raster)

# Define file paths and output names
file_paths <- list(
  bathy = "Variables/pacific/tutuila/Bathymetry_NGDC_AmericanSamoa_90m_all_units.nc",
  sst = "Variables/pacific/tutuila/Tutuila_Sea_Surface_Temperature_jplMUR_Daily_mean_2002-06-02_2023-09-30.nc",
  chla = "Variables/pacific/tutuila/Tutuila_Chlorophyll_A_NPP_VIIRS_Monthly_mean_2012-01-16_2021-01-16.nc",
  dhw = "Variables/pacific/tutuila/Tutuila_Degree_Heating_Weeks_CRW_Daily_DHW.Np10y_1985-03-25_2023-08-31.nc"
)

output_names <- list(
  bathy = "Variables/pacific/tutuila/Processed_Bathy.nc",
  sst = "Variables/pacific/tutuila/Processed_SST.nc",
  chla = "Variables/pacific/tutuila/Processed_Chla.nc",
  dhw = "Variables/pacific/tutuila/Processed_DHW.nc"
)

# Load and process rasters
bathy <- rast(file_paths$bathy)
bathy[bathy <= -100] <- NA
writeCDF(bathy, "Variables/pacific/tutuila/Processed_Bathy.nc")

for (name in names(file_paths)[-1]) {
  
  raster_path <- file_paths[[name]]
  output_path <- output_names[[name]]
  
  raster <- rast(raster_path)
  raster <- terra::resample(raster, bathy)
  raster <- mask(raster, bathy)
  plot(raster)
  
  # Export as NetCDF
  writeCDF(raster, output_path)
}

