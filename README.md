## Scripts for Creating Tiled *h* and *dh/dt* Products from ICESat-2 Altimetry Data

These are instructions for how to create netCDF formated files of the ICESat-2 gridded height products over land ice:
* ATL14 (Antarctic and Greenland ice-sheet surface height)
* ATL15 (Antarctic and Greenland ice-sheet surface height change)

Files ATL14_write2nc.py and ATL15_write2nc.py are located in the surfaceChange/scripts subdirectory.  
They require an ascii input file, containing directory paths, and parameters for domain (x,t) and 
averaging, etc. as shown [here](https://gist.github.com/suzanne64/9483ec8cb8f77200dac2062b3a6da428).

Command line syntax:

\>> python3 \<pathto\>/ATL14_write2nc.py @\<pathto\>/input_args.txt

### Description of the arguments listed in the input_args.txt

The data files to be converted are in hdf5 format. The heights for the ATL14 product are in z0.h5. 
The surface change values for the ATL15 product are in several files, with name starting with dz and 
can include resolution and time lag. These files are to be found in the path in the -b argument. 



