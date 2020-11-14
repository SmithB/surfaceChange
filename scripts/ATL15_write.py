#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:45:47 2020

@author: ben05
"""
import ATL11
import numpy as np
from scipy import stats
import sys, os, h5py, glob, csv
import io
import pointCollection as pc
from utils import projection_details

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import cartopy.crs as ccrs
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
import osgeo.gdal
import datetime as dt
from ATL11.h5util import create_attribute


def ATL15_write(args):
    pathin = os.path.dirname(args.directory)

    def make_dataset(field,data,field_attrs,file_obj,group_obj,scale_dict,dimScale=False):
        dimensions = field_attrs[field]['dimensions'].split(',')
        if field_attrs[field]['datatype'].startswith('int'):
            data = np.nan_to_num(data,nan=np.iinfo(np.dtype(field_attrs[field]['datatype'])).max)
            fillvalue = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
        elif field_attrs[field]['datatype'].startswith('float'):
            data = np.nan_to_num(data,nan=np.finfo(np.dtype(field_attrs[field]['datatype'])).max)
            fillvalue = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
        dset = group_obj.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
        for ii,dim in enumerate(dimensions):
            dset.dims[ii].label = scale_dict[dim.strip()]
            if dimScale:
                dset.make_scale(field)
            else:
                if dim.strip().startswith('Nt'):
                    dset.dims[ii].attach_scale(file_obj[scale_dict[dim.strip()]])
                else:
                    dset.dims[ii].attach_scale(group_obj[scale_dict[dim.strip()]])

        for attr in attr_names:
             if 'dimensions' not in attr and 'datatype' not in attr:
                 create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
        if field_attrs[field]['datatype'].startswith('int'):
            dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
        elif field_attrs[field]['datatype'].startswith('float'):
            dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
        return file_obj

    dz_dict ={'year':'t',     # ATL15product var : mosaic var          
              'year_lag1':'t',
              'year_lag4':'t',
#              'delta_time':'t',
#              'delta_time_lag1':'t',
#              'delta_time_lag4':'t',
              'x':'x',
              'y':'y',
              'cell_area':'cell_area',
              'ice_mask':'mask',
              'data_count':'count',
              'misfit_rms':'misfit_rms',
              'misfit_scaled_rms':'misfit_scaled_rms',
              'delta_h':'dz',
              'delta_h_sigma':'sigma_dz',
              
              'dhdt_lag1':'dzdt_lag1',
              'dhdt_lag4':'dzdt_lag4',
              'dhdt_lag8':'dzdt_lag8',
#              'dhdt':'avg_dzdt',
              'dhdt_lag1_sigma':'sigma_dzdt_lag1',
              'dhdt_lag4_sigma':'sigma_dzdt_lag4',
              'dhdt_lag8_sigma':'sigma_dzdt_lag8',
#              'dhdt_mission':'dhdt_mission',
#              'dhdt_mission_sigma':'dhdt_mission_sigma',
#              'mask_fraction':'mask_fraction',
              }
    scale = {'Nt':'year',
             'Nt_lag1':'year_lag1',
             'Nt_lag4':'year_lag4',
             'Nt_lag8':'year_lag8',
             'Nx':'x',
             'Ny':'y',
             'Nx_10km':'x',
             'Ny_10km':'y',
             'Nx_20km':'x',
             'Ny_20km':'y',
             'Nx_40km':'x',
             'Ny_40km':'y'             
             }
    lags = {
            'file' : ['FH','FH_lag1','FH_lag4'], #,'FH_lag8'],
            'vari' : ['','_lag1','_lag4'], #,'_lag8']
           }
    avgs = ['','_10km','_20km'] #,'_40km']

    # establish output file
    kk=0
    fileout = args.output
    if os.path.isfile(fileout):
        os.remove(fileout)
    with h5py.File(fileout.encode('ASCII'),'w') as fo:
        fo = projection_details(fo,args.Hemisphere)
        # open data attributes file
        with open('../ATL15_output_attrs.csv','r', encoding='utf-8-sig') as attrfile:
            reader=list(csv.DictReader(attrfile))
    
        attr_names=[x for x in reader[0].keys() if x != 'field' and x != 'group']
        
        for kk,ave in enumerate(avgs):
            group_name=f'height_change{ave}'
            field_names = [row['field'] for row in reader if row['group']==group_name]
            # loop over dz*.h5 files for one ave
            for jj in range(len(lags['file'])):
                filein = pathin+'/dz'+ave+lags['vari'][jj]+'.h5'

                if not os.path.isfile(filein):
                    print('No file:',filein)
                    continue
                else:
                    print('Reading file:',filein) #pathin+'/'+os.path.basename(filein))
                lags['file'][jj] = h5py.File(filein,'r')
                dzg=list(lags['file'][jj].keys())[0]
                
                if kk==0:  #establish variables in ROOT
                    for fieldroot in ['year']: #,'delta_time']:
                        field=fieldroot+lags['vari'][jj]
                        data = np.array(lags['file'][jj][dzg][dz_dict[field]])
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
                        if fieldroot == 'year':
                            make_dataset(field,data,field_attrs,fo,fo,scale,dimScale=True)
                        else:
                            make_dataset(field,data,field_attrs,fo,fo,scale,dimScale=False)
    
                if jj==0:
                    gh = fo.create_group(group_name) #'height_change'+ave)
                    # spatial dimension scales for the gh
                    for field in ['x','y']:
                        data = np.array(lags['file'][jj][dzg][dz_dict[field]])
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if (field in row['field'] and row['group']==group_name)}
                        make_dataset(field,data,field_attrs,fo,gh,scale,dimScale=True)
                    
                    for fld in ['delta_h','delta_h_sigma']:
                        field = fld+ave
                        if fld.endswith('sigma'):
                            data = np.array(lags['file'][jj][dzg]['sigma_'+dzg])
                        else:
                            data = np.array(lags['file'][jj][dzg][dzg])
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
                        make_dataset(field,data,field_attrs,fo,gh,scale,dimScale=False)
                    
                    for field in field_names:
                        if not field.startswith('x') and not field.startswith('y') \
                        and not field.startswith('delta_h') and 'lag' not in field:
                            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
                            if dz_dict.get(field)!=None:
                                data = np.array(lags['file'][jj][dzg][dz_dict[field]+ave])
                            else:
                                # place holder data set for now
                                dimensions = field_attrs[field]['dimensions'].split(',')
                                data = np.ndarray(shape=tuple([ii+1 for ii in range(len(dimensions))]),dtype=field_attrs[field]['datatype'])

                            make_dataset(field,data,field_attrs,fo,gh,scale,dimScale=False)
                        
                else:  # one of the lags
                    for fld in ['','_sigma']:
                        field = 'dhdt'+lags['vari'][jj]+fld+ave
                        data = np.array(lags['file'][jj][dzg][dzg])
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
                        make_dataset(field,data,field_attrs,fo,gh,scale,dimScale=False)
                        
            for jj in range(len(lags['file'])):
                lags['file'][jj].close()

    return fileout
    

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--directory','-d', type=str, default=os.getcwd(), help='directory in which to look for mosaicked files')
    parser.add_argument('--output','-o', type=str, default="ATL15.h5", help="output file")
    parser.add_argument('--Hemisphere','-H', type=int, default=1, help='1 for Northern, -1 for Southern')
    args=parser.parse_args()
    print('args',args)
    fileout = ATL15_write(args)



