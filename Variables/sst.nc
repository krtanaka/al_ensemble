CDF       
      time      latitude      	longitude            D   acknowledgement       NOAA Coral Reef Watch Program      cdm_data_type         Grid   comment       �This is a product of the NOAA Coral Reef Watch Daily Global 5km Satellite Coral Bleaching Heat Stress Monitoring Product Suite Version 3.1.    contributor_name      NOAA Coral Reef Watch Program      contributor_role      �Collecting source data and deriving products; performing quality control of products; disseminating, storing, and submitting data to archive   Conventions       CF-1.6, ACDD-1.3, COARDS   creator_email         coralreefwatch@noaa.gov    creator_institution       )NOAA/NESDIS/STAR Coral Reef Watch Program      creator_name      NOAA Coral Reef Watch Program      creator_type      group      creator_url        https://coralreefwatch.noaa.gov/   date_created      2018-03-01T12:00:00Z   date_issued       2022-02-01T21:04:07Z   date_metadata_modified        2020-11-18T12:00:00Z   date_modified         2018-03-01T12:00:00Z   Easternmost_Easting       CJy�   geospatial_bounds         F"POLYGON((-90.0 360.0, 90.0 360.0, 90.0 0.0, -90.0 0.0, -90.0 360.0))"     geospatial_bounds_crs         	EPSG:4326      geospatial_lat_max        A�     geospatial_lat_min        A�33   geospatial_lat_resolution         ?�������   geospatial_lat_units      degrees_north      geospatial_lon_max        CJy�   geospatial_lon_min        CI�    geospatial_lon_resolution         ?�������   geospatial_lon_units      degrees_east   grid_mapping_epsg_code        	EPSG:4326      grid_mapping_inverse_flattening       C� �   grid_mapping_name         latitude_longitude     grid_mapping_semi_major_axis      J¥2   history      	0Tue Jan 23 16:44:33 2024: ncea -v sea_surface_temperature /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202301-0-360.nc /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202302-0-360.nc /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202303-0-360.nc /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202304-0-360.nc /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202305-0-360.nc /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202306-0-360.nc /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202307-0-360.nc /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202308-0-360.nc /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202309-0-360.nc /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202310-0-360.nc /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202311-0-360.nc /mnt/r01/data/CRW/sst/v3.1/monthly/ct5km_sst-mean_v3.1_202312-0-360.nc CRW-sst-2023-clim.nc
Tue Feb  7 07:31:12 2023: ncatted -O -a geospatial_bounds,global,o,c,"POLYGON((-90.0 360.0, 90.0 360.0, 90.0 0.0, -90.0 0.0, -90.0 360.0))" ct5km_sst-mean_v3.1_202301-0-360.nc
Tue Feb  7 07:31:12 2023: ncatted -O -a geospatial_lon_max,global,o,f,359.975 ct5km_sst-mean_v3.1_202301-0-360.nc
Tue Feb  7 07:31:12 2023: ncatted -O -a geospatial_lon_min,global,o,f,0.025 ct5km_sst-mean_v3.1_202301-0-360.nc
Tue Feb  7 07:31:12 2023: ncatted -O -a valid_max,lon,o,f,359.975 ct5km_sst-mean_v3.1_202301-0-360.nc
Tue Feb  7 07:31:12 2023: ncatted -O -a valid_min,lon,o,f,0.025 ct5km_sst-mean_v3.1_202301-0-360.nc
Tue Feb  7 07:31:11 2023: ncap2 -O -s where(lon<0) lon=lon+360 ct5km_sst-mean_v3.1_202301-0-360.nc ct5km_sst-mean_v3.1_202301-0-360.nc
Tue Feb  7 07:31:10 2023: ncks -O --msa_usr_rdr -d lon,0.0,180.0 -d lon,-180.0,0.0 ct5km_sst-mean_v3.1_202301.nc ct5km_sst-mean_v3.1_202301-0-360.nc
This is a product data file of the NOAA Coral Reef Watch Daily Global 5km Satellite Coral Bleaching Heat Stress Monitoring Product Suite Version 3.1 (v3.1) in its NetCDF Version 1.0 (v1.0).
2024-04-10T19:21:12Z (local files)
2024-04-10T19:21:12Z https://oceanwatch.pifsc.noaa.gov/erddap/griddap/CRW_sst_v3_1_2023-clim.nc?sea_surface_temperature%5B(2023-01-31T12:00:00Z)%5D%5B(21.875):(21.025)%5D%5B(201.625):(202.475)%5D&.draw=surface&.vars=longitude%7Clatitude%7Csea_surface_temperature&.colorBar=%7C%7C%7C%7C%7C&.bgColor=0xffffffff   id        /Satellite_Global_5km_SST_Monthly_Mean_Composite    infoUrl       5https://coralreefwatch.noaa.gov/product/5km/index.php      institution       )NOAA/NESDIS/STAR Coral Reef Watch Program      
instrument        �ATSR-1, ATSR-2, AATSR, AVHRR, AVHRR-2, AVHRR-3, VIIRS, GOES Imager, MTSAT Imager, MTSAT 2 Imager, AHI, ABI, SEVIRI, buoy - moored buoy, buoy - drifting buoy, buoy - TAO buoy, surface seawater intake     instrument_vocabulary         *NOAA NODC Ocean Archive System Instruments     keywords     I5km, analysed, composite, coral, coraltemp, crw, data, earth, Earth Science > Oceans > Ocean Temperature > Sea Surface Temperature, Earth Science > Oceans > Ocean Temperature > Water Temperature, Earth Science > Spectral/Engineering > Infrared Wavelengths > Thermal Infrared, engineering, environmental, global, information, infrared, mean, month, monthly, national, nesdis, noaa, ocean, oceans, program, reef, satellite, science, sea, sea_surface_temperature, seawater, service, spectral, spectral/engineering, sst, star, surface, temperature, thermal, time, watch, water, wavelengths      keywords_vocabulary       GCMD Science Keywords      license      KThe data produced by Coral Reef Watch are available for use without restriction, but Coral Reef Watch relies on the ethics and integrity of the user to ensure that the source of the data and products is appropriately cited and credited. When using these data and products, credit and courtesy should be given to NOAA Coral Reef Watch. Please include the appropriate DOI associated with this dataset in the citation. For more information, visit the NOAA Coral Reef Watch website: https://coralreefwatch.noaa.gov. Recommendations for citing and providing credit are provided at https://coralreefwatch.noaa.gov/satellite/docs/recommendations_crw_citation.php. Users are referred to the footer section of the Coral Reef Watch website (https://coralreefwatch.noaa.gov/index.php) for disclaimers, policies, notices pertaining to the use of the data.    metadata_link         5https://coralreefwatch.noaa.gov/product/5km/index.php      naming_authority      gov.noaa.coralreefwatch    NCO       `netCDF Operators version 4.7.5 (Homepage = http://nco.sf.net, Code = https://github.com/nco/nco)   nco_openmp_thread_number            Northernmost_Northing         A�     platform     UShips, drifting buoys, moored buoys, TOGA-TAO buoy arrays, GOES-8 satellite, GOES-9 satellite, GOES-10 satellite, GOES-11 satellite, GOES-12 satellite, GOES-13 satellite, GOES-14 satellite, GOES-15 satellite, GOES-16 satellite, MTSAT-1R satellite, MTSAT-2 satellite, Himawari-8 satellite, Meteosat-8 satellite, Meteosat-9 satellite, Meteoset-10 satellite, Meteosat-11 satellite, Suomi NPP, MetOp-A satellite, MetOp-B satellite, NOAA-9 satellite, NOAA-11 satellite, NOAA-12 satellite, NOAA-14 satellite, NOAA-15 satellite, NOAA-16 satellite, NOAA-17 satellite, NOAA-18 satellite, NOAA-19 satellite.      platform_vocabulary       (NOAA NODC Ocean Archive System Platforms   processing_level      9Derived from L4 satellite sea surface temperaure analysis      product_version       3.1    program       NOAA Coral Reef Watch Program      project       NOAA Coral Reef Watch Program      publisher_email       coralreefwatch@noaa.gov    publisher_institution         )NOAA/NESDIS/STAR Coral Reef Watch Program      publisher_name        NOAA Coral Reef Watch Program      publisher_type        group      publisher_url         https://coralreefwatch.noaa.gov    
references        5https://coralreefwatch.noaa.gov/product/5km/index.php      source        YCoral Reef Watch Daily Global 5km Satellite Sea Surface Temperature v3.1 (CoralTemp v3.1)      	sourceUrl         (local files)      Southernmost_Northing         A�33   spatial_resolution        0.05 degrees   standard_name_vocabulary      CF Standard Name Table v27     summary       �NOAA Coral Reef Watch Global 5km Satellite Sea Surface Temperature (CoralTemp) Monthly Mean Composite for January 2022. This is a product of NOAA Coral Reef Watch Daily Global 5km Satellite Coral Bleaching Heat Stress Monitoring Product Suite     time_coverage_duration        P1M    time_coverage_end         2023-01-31T12:00:00Z   time_coverage_resolution      P1M    time_coverage_start       2023-01-31T12:00:00Z   title         RSea Surface Temperature, Coral Reef Watch, CoralTemp, v3.1 - Cumulative Mean, 2023     Westernmost_Easting       CI�          time             	   _CoordinateAxisType       Time   actual_range      A��@�   A��@�      axis      T      coverage_content_type         
coordinate     ioos_category         Time   	long_name         +reference time of the last day of the month    standard_name         time   time_origin       01-JAN-1970 00:00:00   units         "seconds since 1970-01-01T00:00:00Z          (�   latitude               _CoordinateAxisType       Lat    actual_range      A�33A�     axis      Y      comment       +equirectangular projection and grid centers    coverage_content_type         
coordinate     ioos_category         Location   	long_name         Latitude   standard_name         latitude   units         degrees_north      	valid_max         B��3   	valid_min         ³�3      H  (�   	longitude                  _CoordinateAxisType       Lon    actual_range      CI� CJy�   axis      X      comment       +equirectangular projection and grid centers    coverage_content_type         
coordinate     ioos_category         Location   	long_name         	Longitude      standard_name         	longitude      units         degrees_east   	valid_max         C���   	valid_min         <���      H  )$   sea_surface_temperature                    
   
_FillValue        �         colorBarMaximum       @@         colorBarMinimum                  coverage_content_type         physicalMeasurement    ioos_category         Temperature    	long_name          analysed sea surface temperature   standard_name         sea_surface_temperature    units         degree_C   	valid_max         @��        	valid_min         �i           
   )lA��@�   A�  A���A�33A���A�ffA�  A���A�33A���A�ffA�  A���A�33A���A�ffA�  A���A�33CI� CI��CI��CI�fCI�3CI� CI��CI��CJfCJ3CJ  CJ,�CJ9�CJFfCJS3CJ` CJl�CJy�@9�N���@9��Q�@9�%��X�@9��6�i@9��Q�@9��t�A@9}:Ӡm@9zt�@�@9x�����@9vӠm:@9t�N��@9si�6�@9q��Q@9pm:ӡ@9n���O@9l_���,@9jӠm9@9hQ��@9��Q�@9�K�@9�wwwwx@9�@9��~L@9��N��@9�    @9|(�]@9yb��/�@9w
=p��@9u�,_��@9s�m:�@9rX�%��@9qG�z�@9o\(�@9m�@�t@9k��Q�@9j=p��@9���
=q@9��\(@9�@9�,_���@9�G�z�@9�b��/�@9�""""!@9|_���,@9x�%��Y@9vffffg@9u�Q�@9tDDDDD@9si�6�@9q��S@9p6�i�@9n�Q�@9mi�6�@9k��Q�@9��/�b�@9�z�G�@9�N���@9�""""!@9�N���@9������@9��@�t@9{�����@9w
=p��@9s�
=p�@9s33333@9s33334@9r��/�d@9rX�%��@9q��N�@9pm:ӡ@9n�Q�@9m�@�t@9ȿ%��Y@9����,a@9��Q�@9��6�k@9�K�@9���,_�@9�@�t�@9z�G�|@9s33334@9p��
=p@9p��
=p@9q��Q@9r�\(��@9r��,_�@9r��,_�@9q��N�@9qG�z�@9qG�z�@9��@�u@9�,_���@9������@9�\(�@9�G�z�@9��i�8@9�=p���      �      @9m�@�t@9nz�G�@9o\(�@9qG�z�@9r""""#@9si�6�@9t�@�t@9t�~K�@9u�Q�@9�\(�@9�z�G�@9�ffffg@9�\(�@9�33333@9������@9�������      �      @9k��Q�@9k��Q�@9mi�6�@9n�����@9qG�z�@9s�
=p�@9v/�b��@9xN���@9y�����@9�,_���@9�m:ӡ@9���,_��      �      �      �      �      �      �      @9k�����@9l(�]@9n���O@9p�t�A@9u��X�%@9y,_���@9}i�6�@9���,`@9�i�6�@9�i�6�@9��
=p�@9�
=p���      �      �      �      �      �      �      @9nK�|@9p6�i�@9r�\(��@9wwwwww@9}i�6�@9��\(��@9�ffffg@9�(�[@9������@9��,_��@9�i�6�@9��
=p�      �      �      �      �      �      @9p6�i�@9s33334@9u\(�@9{�X�%�@9��\(��@9�,_���@9����M@:�\(��@9�\(�@9�N���@9�Q��@9��~K��      �      �      �      �      �      �      �      @9|(�]@9���Q@9������@9�~K�@9�@�t�@:�z�H@:\(�@:G�z�@9�K�@9�t�@�@9��~K��      �      �      �      �      �      �      @9�,_���@9�p��
?@9�DDDDD@9������@9�@:
�G�{@:
Ӡm9@:�i�7@:�\(��@9������@9�,_���@9��
=p�@9��Q�@9��N�@9�UUUUU@9�Ӡm:@9������@9��t�@@9��~K@9������@9��t�A@9�N���@9��/�b�@:�����@:p��
=@:
�����@:/�b��@:""""!@9��/�b�@9�/�b��@9�Q� @9�UUUUW@9�Ӡm<@9�i�6�@9�:Ӡm@9�""""#@9�:Ӡm@9�:Ӡm@9��t�A@9�\(�@9������@:�t�A@:\(�@::Ӡm@:	�����@:�~K�@9�\(�@9��%��Y@9��t�@@9�z�H@9�i�6�@9�~K�@9Ȉ����@9�@9����M@9������@9�G�z�@9àm:�@9�z�G�@:��N�@:G�z�@:���P@:(�]@:@�t�@:"""" @9�(�\@9��Q�@9������@9�m:�@9ۅ�Q�@9�UUUUT@9У�
=p@9������@9θQ� @9ϒ��,_@9�m:ӡ@9�6�i�@:33331@:""""!@:     @:�����@:�\)@:�~K�@9��Q�@9��%��Y@9�@9�Ӡm;@9��
=p�@9�     @9�:Ӡm@9�_���,@9ۅ�Q�@9��X�%�@9ۻ����@9��~L@:i�6�@:��N�@:     @:�@�t@:	�6�k@:/�b��@:@9������@9�UUUUW@9�%��X�@9�N���@9興���@9�
=p��@9��,_��@9�UUUUW@9��N��@9��N��@9�~K�