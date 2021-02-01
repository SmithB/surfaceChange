#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:45:47 2020

@author: ben05
"""
import numpy as np
import  os, h5py,  csv
import importlib.resources
from netCDF4 import Dataset

def ATL15_write2nc(args):

    def make_dataset(field,data,field_attrs,file_obj,group_obj,scale,nctype,dimScale=False):
        dimensions = field_attrs[field]['dimensions'].split(',')
        dimensions = tuple(x.strip() for x in dimensions)
        if field_attrs[field]['datatype'].startswith('int'):
            fill_value = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
        elif field_attrs[field]['datatype'].startswith('float'):
            fill_value = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
        data = np.nan_to_num(data,nan=fill_value)

        if dimScale:
            group_obj.createDimension(field_attrs[field]['dimensions'],data.shape[0])

        dsetvar = group_obj.createVariable(field,
                                           nctype[field_attrs[field]['datatype']],
                                           dimensions,
                                           fill_value=fill_value)
            
        dsetvar[:] = data
        for attr in attr_names:
            dsetvar.setncattr(attr,field_attrs[field][attr])
        # add attribute for projection
        if not field.startswith('time'):
            dsetvar.setncattr('grid_mapping','polar_projection')

        return file_obj
    
    dz_dict ={'time':'t',     # for non-lagged vars. {ATL15 outgoing var name: hdf5 incoming var name}          
              'time_lag1':'t',
              'time_lag4':'t',
              'time_lag8':'t',
              'x':'x',
              'y':'y',
              'cell_area':'cell_area',
              'ice_mask':'mask',
              'data_count':'count',
              'misfit_rms':'misfit_rms',
              'misfit_scaled_rms':'misfit_scaled_rms',
              'delta_h':'dz',
              'delta_h_sigma':'sigma_dz',
              'delta_h_10km':'avg_dz_10000m',
              'delta_h_sigma_10km':'sigma_avg_dz_10000m',
              'delta_h_20km':'avg_dz_20000m',
              'delta_h_sigma_20km':'sigma_avg_dz_20000m',
              'delta_h_40km':'avg_dz_40000m',
              'delta_h_sigma_40km':'sigma_avg_dz_40000m',
              }
    nctype = {'float64':'f8',
              'float32':'f4',
              'int8':'i1'}
    scale = {'Nt':'time',
             'Nt_lag1':'time_lag1',
             'Nt_lag4':'time_lag4',
             'Nt_lag8':'time_lag8',
             'Nx':'x',
             'Ny':'y',
             }
    for avg in ['_10km', '_20km', '_40km']:
        for dim in ['x','y']:
            scale[f'N{dim}{avg}']=f'height_change{avg}/{dim}'

    lags = {
            'file' : ['FH','FH_lag1','FH_lag4','FH_lag8'],
            'vari' : ['','_lag1','_lag4','_lag8']
           }
    avgs = ['','_10km','_20km','_40km']

    # establish output file
    fileout = args.base_dir.rstrip('/') + '/ATL15_' + args.region + '_' + args.cycles + '_' + args.Release + '_' + args.version +'.nc'
    print('output file:',fileout)

    with Dataset(fileout,'w',clobber=True) as nc:
        # open data attributes file
        with importlib.resources.path('surfaceChange','resources') as pp:
            with open(os.path.join(pp,'ATL15_output_attrs.csv'),'r', encoding='utf-8-sig') as attrfile:
                reader=list(csv.DictReader(attrfile))

        attr_names=[x for x in reader[0].keys() if x != 'field' and x != 'group']

        for kk,ave in enumerate(avgs):
            # loop over dz*.h5 files for one ave
            for jj in range(len(lags['file'])):
                filein = args.base_dir.rstrip('/')+'/dz'+ave+lags['vari'][jj]+'.h5'
                if not os.path.isfile(filein):
                    print('No file:',args.base_dir.rstrip('/')+'/'+os.path.basename(filein))
                    continue
                else:
                    print('Reading file:',args.base_dir.rstrip('/')+'/'+os.path.basename(filein))
                lags['file'][jj] = h5py.File(filein,'r')  # file object
                dzg=list(lags['file'][jj].keys())[0]      # dzg is group in input file
    
                if jj==0:  # no lag
                    nc.createGroup('height_change'+ave)
                    # each group needs a projection variable
                    if args.region in ['AK','CN','CS','GL','IC','SV','RU']:
                        proj_var = nc.groups['height_change'+ave].createVariable('polar_projection',np.int32,())
                        proj_var.grid_mapping_name = 'Polar Stereographic'
                        proj_var.straight_vertical_longitude_from_pole = -45.0
                        proj_var.latitude_of_projection_origin = 90.0
                        proj_var.standard_parallel = 70.0
                        proj_var.crs_wkt = 'PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3413"]]'
                    elif args.region == 'AA':
                        proj_var = nc.groups['height_change'+ave].createVariable('polar_projection',np.int32,())
                        proj_var.grid_mapping_name = 'Polar Stereographic'
                        proj_var.straight_vertical_longitude_from_pole = 0.0
                        proj_var.latitude_of_projection_origin = -90.0
                        proj_var.standard_parallel = -71.0
                        proj_var.crs_wkt = 'PROJCS["WGS 84 / Antarctic Polar Stereographic",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",-71],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","3031"]]'
                                       
                    # dimension scales for each group
                    for field in ['x','y','time']:
                        data = np.array(lags['file'][jj][dzg][dz_dict[field]])
                        if field == 'time':    # convert to decimal days from 1/1/2018
                            data = (data-2018.)*365.25
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                        make_dataset(field,data,field_attrs,nc,nc.groups['height_change'+ave],scale,nctype,dimScale=True)
                    
                    for fld in ['cell_area','delta_h','delta_h_sigma','ice_mask','data_count','misfit_rms','misfit_scaled_rms']:  # fields that can be ave'd but not lagged
                        if kk>0 and (fld.startswith('misfit') or fld=='ice_mask' or fld=='data_count'): # not in ave'd groups
                            break
                        field = fld+ave
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                        if fld.startswith('delta_h'):  # fields with complicated name changes
                            data = np.array(lags['file'][jj][dzg][dz_dict[field]])
                        else:
                            data = np.array(lags['file'][jj][dzg][dz_dict[fld]])

                        if fld == 'cell_area':
                            data[data==0.0] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
                        make_dataset(field,data,field_attrs,nc,nc.groups['height_change'+ave],scale,nctype,dimScale=False)

                else:  # one of the lags
                    field = 'time'+lags['vari'][jj]
                    data = np.array(lags['file'][jj][dzg]['t'])
                    # convert to decimal days from 1/1/2018
                    data = (data-2018.)*365.25
                    field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                    make_dataset(field,data,field_attrs,nc,nc.groups['height_change'+ave],scale,nctype,dimScale=True)
                    
                    field = 'dhdt'+lags['vari'][jj]+ave
                    data = np.array(lags['file'][jj][dzg][dzg])
                    field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                    make_dataset(field,data,field_attrs,nc,nc.groups['height_change'+ave],scale,nctype,dimScale=False)
                    
                    field = 'dhdt'+lags['vari'][jj]+'_sigma'+ave
                    data = np.array(lags['file'][jj][dzg]['sigma_'+dzg])
                    field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                    make_dataset(field,data,field_attrs,nc,nc.groups['height_change'+ave],scale,nctype,dimScale=False)
                        
            for jj in range(len(lags['file'])):
                try:
                    lags['file'][jj].close()
                except:
                    pass

    return fileout
    

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
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
    args=parser.parse_args()
    print('args',args)
    fileout = ATL15_write2nc(args)



