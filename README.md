## Scripts for Creating Tiled *h* and *dh/dt* Products from ICESat-2 Altimetry Data

These are instructions for how to create netCDF formated files of the ICESat-2 gridded height products
over land ice:
* ATL14 (Antarctic and Greenland ice-sheet surface height)
* ATL15 (Antarctic and Greenland ice-sheet surface height change)

Files ATL14_write2nc.py and ATL15_write2nc.py are located in the surfaceChange/scripts subdirectory.  
They require an ascii input file, containing directory paths, and parameters for domain (x,t) and 
averaging, etc., [input_args.txt](https://gist.github.com/suzanne64/9483ec8cb8f77200dac2062b3a6da428).

Command line syntax:

\>> python3 \<pathto\>/ATL14_write2nc.py @\<pathto\>/input_args.txt
\>> python3 \<pathto\>/ATL15_write2nc.py @\<pathto\>/input_args.txt

### Description of the arguments listed in the input_args.txt that are used for writing netCDF files

-b        base path to height and surface change data files. These are the hdf5 files to be converted netCDF. This directory should contain z0.h5 for ATL14 product. The surface change values for the ATL15 product are in several files, with names starting with dz and which can include resolution and time lag in the filename, ie. dz_10km_lag1.h5. This is also where the output file will be written.


-tiles    path to center tile files

-list11   path to ATL11 files used to make the data files in -b. The ATL11 filenames are included in the metadata.

--region  two character abbreviation for region name. Used to determine gridding projection parameters and part of output file name

--cycles  four digit integer indicating beginning cycle and ending cycle

--Release  three digit integer indicating release of ATL11 data 

--version  two digit integer indicating version of software (?)

### Output file formats

ATL14_rr_c1c2_100m_rel_vv.nc    where rr is the region, c1 is beginning cycle, c2 is ending cycle, rel is release and vv is version
ATL15_rr_c1c2_xxkm_rel_vv.nc    where xx indicates grid spacing in km. Currently there are four output files including 1, 10, 20, 40km grid spacing.

### Size estimates of output files

|   | Greenland | Svalbard |
|----|-----------|----------|
|ATL14_rr_c1c2_100m_rel_vv.nc | 1050 Mb   |  28 Mb|
|ATL15_rr_c1c2_01km_rel_vv.nc |    |  9.2 Mb|
|ATL15_rr_c1c2_10km_rel_vv.nc |    |  0.31 Mb|
|ATL15_rr_c1c2_20km_rel_vv.nc |    |  0.23 Mb|
|ATL15_rr_c1c2_40km_rel_vv.nc |    |  0.20 Mb|



