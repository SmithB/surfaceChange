## Scripts for Creating Tiled *h* and *dh/dt* Products from ICESat-2 Altimetry Data

These are instructions for how to create netCDF formated files of the ICESat-2 gridded height products over land ice:
* ATL14 (Antarctic and Greenland ice-sheet surface height)
* ATL15 (Antarctic and Greenland ice-sheet surface height change)

The scripts are called ATL14_write2nc.py and ATL15_write2nc.py.  They require an arguments input file, with the name inputs_

To run An [example](https://gist.github.com/suzanne64/9483ec8cb8f77200dac2062b3a6da428.js) is here.

Command line syntax:
\>> python3 pathto/ATL14_write2nc.py @pathto/input_args_rr.txt

### Description of the arguements listed in the input_args_rr.txt

