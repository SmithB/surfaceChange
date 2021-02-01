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
        for field in ['x', 'y']: 
            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
            dimensions = field_attrs[field]['dimensions'].split(',')
            fill_value = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
            data = np.array(FH['z0'][dz_dict[field]])
            
            data = np.nan_to_num(data,nan=fill_value)

            nc.createDimension(field_attrs[field]['dimensions'],data.shape[0])
            dsetvar = nc.createVariable(field,
                                        nctype[field_attrs[field]['datatype']],
                                        (field_attrs[field]['dimensions'],),
                                        fill_value=fill_value)
            dsetvar[:] = data
            for attr in attr_names:
                dsetvar.setncattr(attr,field_attrs[field][attr])
        
        for field in [item for item in field_names if item != 'x' and item != 'y']:
            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
            dimensions = field_attrs[field]['dimensions'].split(',')
            dimensions = tuple(x.strip() for x in dimensions)
            data = np.array(FH['z0'][dz_dict[field]])
            if field_attrs[field]['datatype'].startswith('int'):
                fill_value = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
            elif field_attrs[field]['datatype'].startswith('float'):
                fill_value = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
            data = np.nan_to_num(data,nan=fill_value)
            dsetvar = nc.createVariable(field,
                                        nctype[field_attrs[field]['datatype']],
                                        dimensions,
                                        fill_value=fill_value)
            dsetvar[:] = data
            for attr in attr_names:
                dsetvar.setncattr(attr,field_attrs[field][attr])
                
        FH.close()
   
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
    fileout = ATL14_write2nc(args)



