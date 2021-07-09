#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:45:47 2020

@author: ben05
"""

import numpy as np
from scipy import stats
import sys, os, h5py, glob, csv
import io
import pointCollection as pc
import importlib.resources
from netCDF4 import Dataset

from ATL11.h5util import create_attribute

def ATL14_write2nc(args):    
    dz_dict ={'x':'x',   # ATL14 varname : z0.h5 varname
              'y':'y',
              'h':'z0',
              'h_sigma':'sigma_z0',
              'cell_area':'cell_area',
              'ice_mask':'mask',
              'data_count':'count',
              'misfit_rms':'misfit_rms',
              'misfit_scaled_rms':'misfit_scaled_rms',
              }
    nctype = {'float64':'f8',
              'float32':'f4',
              'int8':'i1'}

    # establish output file
    fileout = args.base_dir.rstrip('/') + '/ATL14_' + args.region + '_' + args.cycles + '_' + args.Release + '_' + args.version +'.nc'
    print('output file:',fileout)
   
    with Dataset(fileout,'w',clobber=True) as nc:
        nc.setncattr('GDAL_AREA_OR_POINT','Area')
        nc.setncattr('Conventions','CF-1.6')

        if args.region in ['AK','CN','CS','GL','IC','SV','RU']:
            crs_var = nc.createVariable('Polar_Stereographic',np.byte,())
            crs_var.standard_name = 'Polar_Stereographic'
            crs_var.grid_mapping_name = 'polar_stereographic'
            crs_var.straight_vertical_longitude_from_pole = -45.0
            crs_var.latitude_of_projection_origin = 90.0
            crs_var.standard_parallel = 70.0
            crs_var.scale_factor_at_projection_origin = 1.
            crs_var.false_easting = 0.0
            crs_var.false_northing = 0.0
            crs_var.semi_major_axis = 6378.137
            crs_var.semi_minor_axis = 6356.752
            crs_var.inverse_flattening = 298.257223563
            crs_var.spatial_epsg = '3413'
            crs_var.spatial_ref = 'PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3413"]]'
            crs_var.crs_wkt = ('PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3413"]]')
        elif args.region == 'AA':
            crs_var = nc.createVariable('Polar_Stereographic',np.byte,())
            crs_var.standard_name = 'Polar_Stereographic'
            crs_var.grid_mapping_name = 'polar_stereographic'
            crs_var.straight_vertical_longitude_from_pole = 0.0
            crs_var.latitude_of_projection_origin = -90.0
            crs_var.standard_parallel = -71.0
            crs_var.scale_factor_at_projection_origin = 1.
            crs_var.false_easting = 0.0
            crs_var.false_northing = 0.0
            crs_var.semi_major_axis = 6378.137
            crs_var.semi_minor_axis = 6356.752
            crs_var.inverse_flattening = 298.257223563
            crs_var.spatial_epsg = '3031'
            crs_var.spatial_ref = 'PROJCS["WGS 84 / Antarctic Polar Stereographic",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",-71],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","3031"]]'
            crs_var.crs_wkt = ('PROJCS["WGS 84 / Antarctic Polar Stereographic",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",-71],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","3031"]]')

        # get handle for input file with ROOT and height_change variables.
        FH = h5py.File(args.base_dir.rstrip('/')+'/z0.h5','r')
        if 'z0' not in FH:
            print('no z0.h5 file')
            FH.close()
            exit(-1)
        else:
            print('Reading file:',args.base_dir.rstrip('/')+'/z0.h5')
            
        with importlib.resources.path('surfaceChange','resources') as pp:
            with open(os.path.join(pp,'ATL14_output_attrs.csv'),'r', encoding='utf-8-sig') as attrfile:
                reader=list(csv.DictReader(attrfile))
    
        attr_names=[x for x in reader[0].keys() if x != 'field' and x != 'group']

        field_names = [row['field'] for row in reader if 'ROOT' in row['group']]

        # create dimensions
        for field in ['x', 'y', 'cell_area']: 
            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
            dimensions = field_attrs[field]['dimensions'].split(',')
            dimensions = tuple(x.strip() for x in dimensions)
            fill_value = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
            data = np.array(FH['z0'][dz_dict[field]])
            if field == 'cell_area':
                data[data==0.0] = np.nan 
                cell_area_mask = data  # where cell_area is invalid, so are h and h_sigma
            if field == 'x':
                nc.createDimension(field_attrs[field]['dimensions'],data.shape[0])
                x = data
                xll = np.min(x)
                dx = x[1]-x[0]
            if field == 'y':
                nc.createDimension(field_attrs[field]['dimensions'],data.shape[0])
                y = data
                yll = np.max(y)
                dy = y[0]-y[1]    
                
            data = np.nan_to_num(data,nan=fill_value)
            dsetvar = nc.createVariable(field,
                                        nctype[field_attrs[field]['datatype']],
                                        dimensions, zlib=True, least_significant_digit=4,
                                        fill_value=fill_value)
            dsetvar[:] = data
            for attr in attr_names:
                dsetvar.setncattr(attr,field_attrs[field][attr])
            # add attributes for projection
            if field == 'x':
                dsetvar.standard_name = 'projection_x_coordinate'
            if field == 'y':
                dsetvar.standard_name = 'projection_y_coordinate'                
        crs_var.GeoTransform = (xll,dx,0,yll,0,dy)
                
        for field in [item for item in field_names if item != 'x' and item != 'y' and item != 'cell_area']:
            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
            dimensions = field_attrs[field]['dimensions'].split(',')
            dimensions = tuple(x.strip() for x in dimensions)
            data = np.array(FH['z0'][dz_dict[field]])
            if field.startswith('h'):
                data[np.isnan(cell_area_mask)] = np.nan
            if field_attrs[field]['datatype'].startswith('int'):
                fill_value = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
            elif field_attrs[field]['datatype'].startswith('float'):
                fill_value = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
            data = np.nan_to_num(data,nan=fill_value)
            dsetvar = nc.iable(field,
                                        nctype[field_attrs[field]['datatype']],
                                        dimensions, zlib=True, least_significant_digit=4,
                                        fill_value=fill_value)
            dsetvar[:] = data
            for attr in attr_names:
                dsetvar.setncattr(attr,field_attrs[field][attr])
            dsetvar.setncattr('grid_mapping','Polar_Stereographic')
                
        FH.close()
   
    return fileout
    
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, fromfile_prefix_chars='@')
    parser.add_argument('-b','--base_dir', type=str, default=os.getcwd(), help='directory in which to look for mosaicked .h5 files')
    parser.add_argument('-rr','--region', type=str, help='2-letter region indicator \n'
                                                         '\t AA: Antarctica \n'
                                                         '\t AK: Alaska \n'
                                                         '\t CN: Arctic Canada North \n'
                                                         '\t CS: Arctic Canada South \n'
                                                         '\t GL: Greeland and peripheral ice caps \n'
                                                         '\t IC: Iceland \n'
                                                         '\t SV: Svalbard \n'
                                                         '\t RU: Russian Arctic')
    parser.add_argument('-c','--cycles', type=str, help="4-digit number specifying first/last cycles for output filename")
    parser.add_argument('-R','--Release', type=str, help="3-digit release number for output filename")
    parser.add_argument('-v','--version', type=str, help="2-digit version number for output filename")
    args, _=parser.parse_known_args()
    print('args',args)
    fileout = ATL14_write2nc(args)



