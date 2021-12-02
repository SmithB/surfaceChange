## Scripts for Creating Tiled *h* and *dh/dt* Products from ICESat-2 Altimetry Data

These are instructions for how to create netCDF formated files of the ICESat-2 gridded height products
over land ice:
* ATL14 (Antarctic and Greenland ice-sheet surface height)
* ATL15 (Antarctic and Greenland ice-sheet surface height change)

Files ATL14_write2nc.py and ATL15_write2nc.py are located in the surfaceChange/scripts subdirectory.  
They require an ascii input file, containing directory paths, and parameters for domain (x,t) and
averaging, [input_args.txt](https://gist.github.com/suzanne64/9483ec8cb8f77200dac2062b3a6da428).

Command line syntax:

```bash
>> python3 \<pathto\>/ATL14_write2nc.py @\<pathto\>/input_args.txt   
>> python3 \<pathto\>/ATL15_write2nc.py @\<pathto\>/input_args.txt
```

### Description of the arguments listed in the input_args.txt that are used for writing netCDF files

`-b`        base path to height and surface height change data files. These are the hdf5 files to be consolidated and converted to netCDF. This directory should contain z0.h5 for the ATL14 product. The surface height change values for the ATL15 product are in several files, with names starting with dz and which can include resolution and time lag in the filename, ie. dz_10km_lag1.h5. This is also where the output netCDF files will be written.

`--region`  two character abbreviation for region name. Used to determine gridding projection parameters and is part of output file name.

`--cycles`  four digit integer indicating beginning cycle and ending cycle

`--Release`  three digit integer indicating release of ATL11 data

`--version`  two digit integer indicating version of software (?)

### Output file formats

The output files are in netCDF format. The filenames are

`ATL14_rr_c1c2_100m_rel_vv.nc`    where `rr` is the region, `c1` is beginning cycle, `c2` is ending cycle, `rel` is release and `vv` is version.  
`ATL15_rr_c1c2_xxkm_rel_vv.nc`    where `xx` indicates grid spacing in km. Currently there are four ATL15 output files including 1, 10, 20, 40km grid spacing.

### Size estimates of output files

|   | Greenland | Svalbard |
|----|-----------|----------|
|ATL14_rr_c1c2_100m_rel_vv.nc | 1050 Mb   |  28 Mb|
|ATL15_rr_c1c2_01km_rel_vv.nc |  343 Mb  |  9.2 Mb|
|ATL15_rr_c1c2_10km_rel_vv.nc |    3.5 Mb|  0.31 Mb|
|ATL15_rr_c1c2_20km_rel_vv.nc |    1.4 Mb|  0.23 Mb|
|ATL15_rr_c1c2_40km_rel_vv.nc |    0.8 Mb|  0.20 Mb|


### Contributing Code or Examples

We follow a standard Forking Workflow for code changes and additions.
Submitted code goes through the pull request process for evaluation and comments.

#### General Guidelines

- Make each pull request as small and simple as possible
- Commit messages should be clear and describe the changes
- Larger changes should be broken down into their basic components and integrated separately
- If possible, bug fixes should be their own pull requests with an associated [GitHub issue](https://github.com/SmithB/surfaceChange/issues)
- Write a descriptive pull request message with a clear title

#### Steps to Contribute

1. Fork the repository to your personal GitHub account by clicking the "Fork" button on the project [main page](https://github.com/SmithB/surfaceChange).  This creates your own server-side copy of the repository.
2. Create a work environment to make your changes by cloning to your local system.
3. Add your fork as the `origin` remote and the original project repository as the `upstream` remote.  While this step isn't a necessary, it allows you to keep your fork up to date in the future.
4. Create a new branch to do your work.
5. Make your changes on the new branch and add yourself to the list of [project contributors](./CONTRIBUTORS.md).
6. Push your work to GitHub under your fork of the project.
7. Submit a [Pull Request](https://github.com/SmithB/surfaceChange/pulls) from your forked branch to the project repository.
