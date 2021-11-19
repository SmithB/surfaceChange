#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:45:47 2020

@author: ben05
"""

import numpy as np
from scipy import stats
import sys, os, h5py, glob, csv
import io, re
import pointCollection as pc
import pkg_resources
from netCDF4 import Dataset
import matplotlib.pyplot as plt

from ATL14_attrs_meta import write_atl14meta
#from ATL11.h5util import create_attribute

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
    fileout = args.base_dir.rstrip('/') + '/ATL14_' + args.region + '_' + args.cycles + '_100m_' + args.Release + '_' + args.version +'.nc'
    print('output file:',fileout)
   
    with Dataset(fileout,'w',clobber=True) as nc:
        nc.setncattr('GDAL_AREA_OR_POINT','Area')
        nc.setncattr('Conventions','CF-1.6')

        if args.region in ['AK','CN','CS','GL','IS','SV','RA']:
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

        # make tile_stats group (ATBD 4.1.2.1, Table 3)
        tilegrp = nc.createGroup('tile_stats')   
        tile_stats = {'x': { 'data': [], 'description':'tile-center x-coordinate, in projected coordinates', 'mapped':np.array(())},
                      'y': { 'data': [], 'description':'tile-center y-coordinate, in projected coordinates', 'mapped':np.array(())},
                      'N_data': { 'data': [], 'description':'number of data used in fit', 'mapped':np.array(())},
                      'RMS_data': { 'data': [], 'description':'root mean of squared, scaled data misfits', 'mapped':np.array(())},
                      'RMS_bias': { 'data': [], 'description':'root mean of squared, scaled bias values', 'mapped':np.array(())},
                      'N_bias': { 'data': [], 'description':'number of bias values solved for', 'mapped':np.array(())},
                      'RMS_d2z0dx2': { 'data': [], 'description':'root mean square of the constraint equation residuals for the second spatial derivative of z0', 'mapped':np.array(())},
                      'RMS_d2zdt2': { 'data': [], 'description':'root mean square of the constraint equation residuals for the second temporal derivative of dz', 'mapped':np.array(())},
                      'RMS_d2zdx2dt' : { 'data': [], 'description':'root mean square of the constraint equation residuals for the second temporal derivative of dz/dt', 'mapped':np.array(())},
                      'sigma_xx0' : { 'data': [], 'description':'Weighting values for the constraint equations on the second spatial derivatives of the DEM', 'mapped':np.array(())},
                      'sigma_tt' : { 'data': [], 'description':'Weighting values for the constraint equations on the second temporal derivatives of the surface height', 'mapped':np.array(())},
                      'sigma_xxt' : { 'data': [], 'description':'Weighting values for the constraint equations on the second spatial derivatives of the height-change rate', 'mapped':np.array(())}
                      }
        
        # work through the tiles in all three subdirectories
        for sub in ['centers','edges','corners']:
            files = os.listdir(os.path.join(args.base_dir,sub))
            files = [f for f in files if f.endswith('.h5')]
            for file in files:
                try:
                    tile_stats['x']['data'].append(int(re.match(r'^.*E(.*)\_.*$',file).group(1)))
                except Exception as e:
                    print(f"problem with [ {file} ], skipping")
                    continue
                tile_stats['y']['data'].append(int(re.match(r'^.*N(.*)\..*$',file).group(1)))
                with h5py.File(os.path.join(args.base_dir,sub,file),'r') as h5:
                    tile_stats['N_data']['data'].append( np.sum(h5['data']['three_sigma_edit'][:]) )
                    tile_stats['RMS_data']['data'].append( h5['RMS']['data'][()] )  # use () for getting a scalar.
                    tile_stats['RMS_bias']['data'].append( np.sqrt(np.mean((h5['bias']['val'][:]/h5['bias']['expected'][:])**2)) )
                    tile_stats['N_bias']['data'].append( len(h5['bias']['val'][:]) )  #### or all BUT the zeros.
                    tile_stats['RMS_d2z0dx2']['data'].append( h5['RMS']['grad2_z0'][()] )
                    tile_stats['RMS_d2zdt2']['data'].append( h5['RMS']['d2z_dt2'][()] )
                    tile_stats['RMS_d2zdx2dt']['data'].append( h5['RMS']['grad2_dzdt'][()] )
                    tile_stats['sigma_xx0']['data'].append( h5['E_RMS']['d2z0_dx2'][()] )
                    tile_stats['sigma_tt']['data'].append( h5['E_RMS']['d2z_dt2'][()] )
                    tile_stats['sigma_xxt']['data'].append( h5['E_RMS']['d3z_dx2dt'][()] )
                    
        # establish output grids
        for key in tile_stats.keys():
            if key == 'x' or key == 'y' or key == 'N_data':  #integers
                tile_stats[key]['mapped'] = np.zeros( [len(np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+40,40)),
                                                       len(np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+40,40))], 
                                                       dtype=int)
            else:
                tile_stats[key]['mapped'] = np.zeros( [len(np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+40,40)),
                                                       len(np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+40,40))],
                                                       dtype=float)

        # fill grids
        for i, (yt,xt) in enumerate(zip(tile_stats['y']['data'],tile_stats['x']['data'])):
            for key in tile_stats.keys():
                # fact helps convert x,y in km to m
                if key != 'x' and key != 'y':
                    tile_stats[key]['mapped'][int((yt-np.min(tile_stats['y']['data']))/40),int((xt-np.min(tile_stats['x']['data']))/40)] = \
                    tile_stats[key]['data'][i]
                    tile_stats[key]['mapped'] = np.ma.masked_where(tile_stats[key]['mapped'] == 0, tile_stats[key]['mapped'])   

        # make dimensions, fill them as variables
        tilegrp.createDimension('y',len(np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+40,40)))
        y = tilegrp.createVariable('y', np.dtype('int32'), ('y',))
        y[:]=np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+40,40) * 1000 # convert from km to meter
        y.units = 'meter'
        y.description = tile_stats['y']['description']
        y.grid_mapping = 'Polar_Stereographic'
        
        tilegrp.createDimension('x',len(np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+40,40)))
        x = tilegrp.createVariable('x', np.dtype('int32'), ('x',))
        x[:]=np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+40,40) * 1000
        x.units = 'meter'
        x.description = tile_stats['x']['description']
        x.grid_mapping = 'Polar_Stereographic'

        # create tile_stats/ variables in .nc file
        for key in tile_stats.keys():
            if key != 'x' and key != 'y':
                if key == 'N_data': 
                    dsetvar = tilegrp.createVariable(key,np.dtype('int32'),('y','x'),fill_value=np.iinfo(np.dtype('int32')).max, zlib=True)
                else:
                    dsetvar = tilegrp.createVariable(key,np.dtype('float64'),('y','x'),fill_value=np.finfo(np.dtype('float64')).max, zlib=True)
                    dsetvar.least_significant_digit = 4
                dsetvar[:] = tile_stats[key]['mapped'][:]
                dsetvar.setncattr('description',tile_stats[key]['description'])
                dsetvar.setncattr('grid_mapping','Polar_Stereographic')

        # get handle for input file with ROOT and height_change variables.
        FH = h5py.File(args.base_dir.rstrip('/')+'/z0.h5','r')
        if 'z0' not in FH:
            print('no z0.h5 file')
            FH.close()
            exit(-1)
        else:
            print('Reading file:',args.base_dir.rstrip('/')+'/z0.h5')
            
        attrFile = pkg_resources.resource_filename('surfaceChange','resources/ATL14_output_attrs.csv') 
        with open(attrFile,'r', encoding='utf-8-sig') as attrfile:
            reader=list(csv.DictReader(attrfile))
    
        attr_names=[x for x in reader[0].keys() if x != 'field' and x != 'group']

        field_names = [row['field'] for row in reader if 'ROOT' in row['group']]

        # create dimensions
        for field in ['x', 'y', 'ice_mask', 'cell_area']: 
            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
            dimensions = field_attrs[field]['dimensions'].split(',')
            dimensions = tuple(x.strip() for x in dimensions)
            if field != 'ice_mask':
                fill_value = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
            else:
                fill_value = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max                
            data = np.array(FH['z0'][dz_dict[field]])
            if field == 'ice_mask':
                ice_mask = data
            if field == 'cell_area':
                data = data*ice_mask
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
                
        for field in [item for item in field_names if item != 'x' and item != 'y' and item != 'ice_mask' and item != 'cell_area']:
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
            dsetvar = nc.createVariable(field,
                                        nctype[field_attrs[field]['datatype']],
                                        dimensions, zlib=True, least_significant_digit=4,
                                        fill_value=fill_value)
            dsetvar[:] = data
            for attr in attr_names:
                dsetvar.setncattr(attr,field_attrs[field][attr])
            dsetvar.setncattr('grid_mapping','Polar_Stereographic')

        ncTemplate = pkg_resources.resource_filename('surfaceChange','resources/atl14_metadata_template.nc')
        write_atl14meta(nc, fileout, ncTemplate, args)

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
                                                         '\t IS: Iceland \n'
                                                         '\t SV: Svalbard \n'
                                                         '\t RA: Russian Arctic')
    parser.add_argument('-c','--cycles', type=str, help="4-digit number specifying first/last cycles for output filename")
    parser.add_argument('-R','--Release', type=str, help="3-digit release number for output filename")
    parser.add_argument('-v','--version', type=str, help="2-digit version number for output filename")
    parser.add_argument('-list11','--ATL11_lineage_dir', type=str, help='directory in which to look for ATL11 .h5 filenames')
    parser.add_argument('-tiles','--tiles_dir', type=str, help='directory in which to look for tile .h5 files, defaults to [base_dir]/centers')
    parser.add_argument('--ATL11_index', type=str, help='GeoIndex file pointing to ATL11 data')
    args, _=parser.parse_known_args()
    
    if args.ATL11_lineage_dir is None:
        # if ATL11 lineage_dir is not specified, assume that the grandparent of the ATL11_index works
        args.ATL11_lineage_dir=os.path.dirname(os.path.dirname(args.ATL11_index))
    
    if args.tiles_dir is None:
        args.tiles_dir=os.path.join(args.base_dir, 'centers')
    
    print('args:',args)
    
    
    
    fileout = ATL14_write2nc(args)
