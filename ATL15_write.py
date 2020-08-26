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
from PointDatabase.mapData import mapData

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
import imageio
import datetime as dt
from ATL11.h5util import create_attribute


def ATL15_write():
    
#    def make_dataset(field,data,field_attrs,group_obj,file_obj,dimScale=False):
#        dimensions = field_attrs[field]['dimensions'].split(',')
#        data = np.nan_to_num(data,nan=np.finfo(np.dtype(field_attrs[field]['datatype'])).max)
#        fillvalue = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
#        dset = file_obj.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
#        for ii,dim in enumerate(dimensions):
#            dset.dims[ii].label = scale[dim.strip()]
#            if fieldroot == 'year':
#                dset.make_scale(field)
#            else:    
#                dset.dims[ii].attach_scale(fo[scale[dim.strip()]])
#        for attr in attr_names:
#             if 'dimensions' not in attr and 'datatype' not in attr:
#                 create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
#        if field_attrs[field]['datatype'].startswith('int'):
#            dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
#        elif field_attrs[field]['datatype'].startswith('float'):
#            dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
#        print(file_obj)

    dz_dict ={'year':'t',               
              'year_lag1':'t',
              'year_lag4':'t',
              'delta_time':'t',
              'delta_time_lag1':'t',
              'delta_time_lag4':'t',
              'x':'x',
              'y':'y',
#              'cell_area':'area',
              'data_count':'count',
              'misfit_rms':'misfit_rms',
              'misfit_scaled_rms':'misfit_scaled_rms',
              'delta_h':'dz',
              'delta_h_sigma':'sigma_dz',
              'dhdt':'dzdt',
#              'dhdt_lag1_sigma':'sigma_avg_dzdt_lag1',
#              'dhdt_lag4':'dhdt_lag4',
#              'dhdt_lag4_sigma':'dhdt_lag4_sigma',
#              'dhdt_mission':'dhdt_mission',
#              'dhdt_mission_sigma':'dhdt_mission_sigma',
#              'ice_mask':'ice_mask',
#              'mask_fraction':'mask_fraction',
              }
#    dz_10km_dict ={'x':'x',
#                   'y':'y',
#                   'data_count':'count',
#                   'delta_h':'avg_dz_10000m',
#                   'delta_h_sigma':'sigma_avg_dz_10000m',
#              }
    scale = {'Nt':'year',
             'Nt_lag1':'year_lag1',
             'Nt_lag4':'year_lag4',
             'Nx':'x',
             'Ny':'y',
             'Nx_10km':'x_10km',
             'Ny_10km':'y_10km',
             'Nx_20km':'x_20km',
             'Ny_20km':'y_20km',
             'Nx_40km':'x_40km',
             'Ny_40km':'y_40km',             
             }
    lags = {
            'lag': ['','_lag1','_lag4'],
            'fileh' : ['FH','FH_lag1','FH_lag4']
           }

    # establish output file
    fileout = 'ATL15_yyyymmdd.h5'
#    print('output file',fileout)
    if os.path.isfile(fileout):
        os.remove(fileout)
    with h5py.File(fileout.encode('ASCII'),'w') as fo:
        # open input files with ROOT and height_change variables.
#            if 'dz' not in FH:
#                print('no dz group')
#                FH.close()
#                exit(-1)
#        # height_change groups
#        gh=f0.greate_
        with open('ATL15_output_attrs_sd.csv','r', encoding='utf-8-sig') as attrfile:
            reader=list(csv.DictReader(attrfile))
        group_names = set([row['group'] for row in reader])
    
        attr_names=[x for x in reader[0].keys() if x != 'field' and x != 'group']
#        print('line 68',attr_names)
#        print('group names',group_names)
        
        # work ROOT group first
        field_names = [row['field'] for row in reader if 'ROOT' in row['group']]
        # open input files for group height_change
        for jj in range(len(lags['fileh'])):
            lags['fileh'][jj] = h5py.File('dz'+lags['lag'][jj]+'.h5','r')
            #establish variables that are dimension scales first
            for fieldroot in ['year','delta_time']:
                field=fieldroot+lags['lag'][jj]
                print('line 108 field',field,lags['lag'][jj],fieldroot,jj)
                dzg=list(lags['fileh'][jj].keys())[0]
#                if dzg.keys() != None:
                data = np.array(lags['fileh'][jj][dzg][dz_dict[field]])
                print(jj,data)
#                if dz_dict.get(field)!=None:   # if key in dz_dict
                
                field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
#                make_dataset(field,data,field_attrs,[],fo)
#                exit(-1)
                dimensions = field_attrs[field]['dimensions'].split(',')
                data = np.nan_to_num(data,nan=np.finfo(np.dtype(field_attrs[field]['datatype'])).max)
                fillvalue = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
                dset = fo.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
                if fieldroot == 'year':
                    dset.make_scale(field)
                for ii,dim in enumerate(dimensions):
                    dset.dims[ii].label = scale[dim.strip()]
                    if fieldroot == 'year':
                        dset.make_scale(field)
                    else:    
                        dset.dims[ii].attach_scale(fo[scale[dim.strip()]])
                for attr in attr_names:
                     if 'dimensions' not in attr and 'datatype' not in attr:
                         create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
                if field_attrs[field]['datatype'].startswith('int'):
                    dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
                elif field_attrs[field]['datatype'].startswith('float'):
                    dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max

            if jj==0:
                gh = fo.create_group('height_change')
                for field in ['x','y']:
                    data = np.array(lags['fileh'][jj][dzg][dz_dict[field]])
                    print('     jj',jj,field,data.shape)
                    field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
                    dimensions = field_attrs[field]['dimensions'].split(',')
                    data = np.nan_to_num(data,nan=np.finfo(np.dtype(field_attrs[field]['datatype'])).max)
                    fillvalue = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
                    dset = gh.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
                    dset.make_scale(field)
                    for ii,dim in enumerate(dimensions):
                        dset.dims[ii].label = scale[dim.strip()]
                    for attr in attr_names:
                         if 'dimensions' not in attr and 'datatype' not in attr:
                             create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
                    if field_attrs[field]['datatype'].startswith('int'):
                        dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
                    elif field_attrs[field]['datatype'].startswith('float'):
                        dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
                for field in ['data_count','misfit_rms','misfit_scaled_rms','delta_h','delta_h_sigma']:
                    # missign: 'cell_area','dhdt_mission','dhdt_mission_sigma','ice_mask'
                    data = np.array(lags['fileh'][jj][dzg][dz_dict[field]])
                    print('     jj',jj,field,data.shape)
                    field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
                    dimensions = field_attrs[field]['dimensions'].split(',')
                    data = np.nan_to_num(data,nan=np.finfo(np.dtype(field_attrs[field]['datatype'])).max)
                    fillvalue = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
                    dset = gh.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
                        
                    for ii,dim in enumerate(dimensions):
                        dset.dims[ii].label = scale[dim.strip()]
                        print(field,ii,dim)
                        if field == 'x' or field =='y':
                            dset.make_scale(field)
                        else:
                            if dim.strip().startswith('Nt'):
                                dset.dims[ii].attach_scale(fo[scale[dim.strip()]])
                            else:
                                dset.dims[ii].attach_scale(gh[scale[dim.strip()]])
                    for attr in attr_names:
                         if 'dimensions' not in attr and 'datatype' not in attr:
                             create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
                    if field_attrs[field]['datatype'].startswith('int'):
                        dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
                    elif field_attrs[field]['datatype'].startswith('float'):
                        dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max

            else:  # one of the lags
                for fld in ['','_sigma']:
                    field = 'dhdt'+lags['lag'][jj]+fld
                    infield = fld+'_'+dz_dict['dhdt']+lags['lag'][jj]
                    infield = infield[1:]
                    print('line 210',field,infield)
                    data = np.array(lags['fileh'][jj][dzg][infield])
                    print('     jj',jj,field,data.shape)
                    field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
                    dimensions = field_attrs[field]['dimensions'].split(',')
                    data = np.nan_to_num(data,nan=np.finfo(np.dtype(field_attrs[field]['datatype'])).max)
                    fillvalue = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
                    dset = gh.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
                        
                    for ii,dim in enumerate(dimensions):
                        dset.dims[ii].label = scale[dim.strip()]
                        print(field,ii,dim)
                        if field == 'x' or field =='y':
                            dset.make_scale(field)
                        else:
                            if dim.strip().startswith('Nt'):
                                dset.dims[ii].attach_scale(fo[scale[dim.strip()]])
                            else:
                                dset.dims[ii].attach_scale(gh[scale[dim.strip()]])
                    for attr in attr_names:
                         if 'dimensions' not in attr and 'datatype' not in attr:
                             create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
                    if field_attrs[field]['datatype'].startswith('int'):
                        dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
                    elif field_attrs[field]['datatype'].startswith('float'):
                        dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max





        for ii in range(len(lags['fileh'])):
            lags['fileh'][ii].close()
#        
#        for field in [item for item in field_names if not item.startswith('year')]:
#            if 'lag1' in field:   # need to fix when we get delta_time_1lag
#                data = np.array(FH['dz'][dz_dict[field]][1:])
#            elif 'lag4' in field:
#                data = np.array(FH['dz'][dz_dict[field]][4:])
#            else:
#                data = np.array(FH['dz'][dz_dict[field]])
#                
#            # read attrs from .csv
#            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
#            dimensions = field_attrs[field]['dimensions'].split(',')
#            data = np.nan_to_num(data,nan=np.finfo(np.dtype(field_attrs[field]['datatype'])).max)
#            fillvalue = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
#            dset = fo.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
#            for ii,dim in enumerate(dimensions):
#                dset.dims[ii].label = scale[dim.strip()]
#                if dim.strip() == 'Nt':
#                    dset.dims[ii].attach_scale(fo[scale[dim.strip()]])
#            for attr in attr_names:
#                 if 'dimensions' not in attr and 'datatype' not in attr:
#                     create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
#            if field_attrs[field]['datatype'].startswith('int'):
#                dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
#            elif field_attrs[field]['datatype'].startswith('float'):
#                dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
#
#        # four height change groups
#        for gp in ['height_change']: #,'height_change_10km','height_change_20km','height_change_40km']:
#            grpstr=gp[-5:]
#            if grpstr=='hange':
#                grpstr = ''
#            print(gp,grpstr)
#            g=fo.create_group(gp)
#            field_names = [row['field'] for row in reader if row['group'] == gp]            
#            print('group', gp,' has these fields',field_names)
#            # open three input files per group
#            
#            FH = h5py.File('dz'+grpstr+'.h5','r')
#            if len(grpstr)>0:
#                FHlag1 = h5py.File('dz'+grpstr+'_lag1.h5','r')
#                FHlag4 = h5py.File('dz'+grpstr+'_lag4.h5','r')
#            # establish dimension scale variables
#            for fld in ['x','y']:
#                field = fld+grpstr
#                print('field',field)
#                data = np.array(FH['dz'][dz_dict[field]])
#                field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
#                dimensions = field_attrs[field]['dimensions'].split(',')
#                data = np.nan_to_num(data,nan=np.finfo(np.dtype(field_attrs[field]['datatype'])).max)
#                fillvalue = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
#                dset = g.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
#                dset.make_scale(field)
#                for ii,dim in enumerate(dimensions):
#                    print(field,dim)
#                    dset.dims[ii].label = scale[dim.strip()]
#                for attr in attr_names:
#                     if 'dimensions' not in attr and 'datatype' not in attr:
#                         create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
#                if field_attrs[field]['datatype'].startswith('int'):
#                    dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
#                elif field_attrs[field]['datatype'].startswith('float'):
#                    dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
#                print()
#
#            for fld in [item for item in field_names if not item.startswith(('x', 'y'))]:
#                field = fld+grpstr
#                # read attrs from .csv
#                field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
#                dimensions = field_attrs[field]['dimensions'].split(',')
#                if dz_dict.get(field)!=None:   # if key in dz_dict
#                    if 'lag1' in field:   
#                        data = np.array(FH['dz'][dz_dict[field]][1:])
#                    elif 'lag4' in field:
#                        data = np.array(FH['dz'][dz_dict[field]][4:])
#                    else:
#                        data = np.array(FH['dz'][dz_dict[field]])
#                else:
#                    data = np.ndarray(shape=tuple([ii+1 for ii in range(len(dimensions))]),dtype=float)
#                    
#                if field_attrs[field]['datatype'].startswith('int'):
#                    data = np.nan_to_num(data,nan=np.iinfo(np.dtype(field_attrs[field]['datatype'])).max)
#                    data = data.astype('int')  # don't change to int before substituting nans with invalid.
#                    fillvalue = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
#                elif field_attrs[field]['datatype'].startswith('float'):
#                    data = np.nan_to_num(data,nan=np.finfo(np.dtype(field_attrs[field]['datatype'])).max)
#                    fillvalue = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
#                dset = g.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
#                for ii,dim in enumerate(dimensions):
#                    dset.dims[ii].label = scale[dim.strip()]
#                    if dim.strip() == 'Nt':
#                        dset.dims[ii].attach_scale(fo[scale[dim.strip()]])
#                for attr in attr_names:
#                     if 'dimensions' not in attr and 'datatype' not in attr:
#                         create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
#                if field_attrs[field]['datatype'].startswith('int'):
#                    dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
#                elif field_attrs[field]['datatype'].startswith('float'):
#                    dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
            
            
            
            
            
            
#                for field in field_names:
#                    field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
#                    dimensions = field_attrs[field]['dimensions'].split(',')
#                    if dz_dict.get(field)!=None:   # if key in dz_dict
#                        if field.endswith('_lag1'):   # need to fix when we get delta_time_1lag
#                            data = np.array(FH['dz'][dz_dict[field]][1:])
#                        elif field.endswith('_lag4'):
#                            data = np.array(FH['dz'][dz_dict[field]][4:])
#                        else:
#                            data = np.array(FH['dz'][dz_dict[field]])
#                    else:
#                        data = np.ndarray(shape=tuple([ii+1 for ii in range(len(dimensions))]),dtype=float)
#                        
#                    dset = g.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
#                    print('dimensions',field,dimensions)
#                    for ii,dim in enumerate(dimensions):
#                        dset.dims[ii].label = dimensions[ii]
#                        if dim is 'Nt':
#                            print('line 144',dim)
#                            print(scale[dim.strip()])
#                            print()
#                            dset.dims[ii].attach_scale(g[scale[dim.strip()]])
                     
                     
#                        .attach_scale(g['cycle_number'])
#                        if dim is 'Nt': 
#                            print(dim)
#                            dset.dims[ii].attach_scale(g['year'])
#                        if dim is 'Nt': 
#                            print(dim)
#                            dset.dims[ii].attach_scale(g['year'])
#                        if dim is 'Nt': 
#                            print(dim)
#                            dset.dims[ii].attach_scale(g['year'])
#                    for attr in attr_names:
#                         if 'dimensions' not in attr and 'datatype' not in attr:
#                             create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
#                    if field_attrs[field]['datatype'].startswith('int'):
#                        dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
#                    elif field_attrs[field]['datatype'].startswith('Float'):
#                        dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
                    
                    
                #dset.dims[0].label = field
                
    #            g.create_dataset('year_lag1',data=np.array(FH['dz']['t'][1:]))
    #            field_attrs = {'year_lag1': {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader}
    #            g.create_dataset('year_lag4',data=np.array(FH['dz']['t'][4:]))
    #            field_attrs = {'year_lag4': {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader}
                
    #            print(dz_dict['t'],np.array(FH['dz']['t']))
    #            g.create_dataset(dz_dict['t'],data=np.array(FH['dz']['t']))
    #            for key,val in FH['dz'].items():   
    #                print('line 71',key,dz_dict[key])
    #                field_attrs = {row[dz_dict[key]]: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader}
    #                print('line 86',field_attrs)
    #                if key is 't':
    #                    g.create_dataset(dz_dict[key],data=np.array(FH['dz'][key]))
    #                else:
    #                    gh.create_dataset(dz_dict[key],data=np.array(FH['dz'][key]))
                        
                        
#        FH.close()
#        FH10.close()
#        FH20.close()
#        FH40.close()
            
            
            
            
            
#            for key,val in attrfile something    FH[pt].items():
#                print('items in ROOT',key,val)
#            dset = g.create_dataset(field.encode('ASCII'),data=data,chunks=True,compression=6,dtype=field_attrs[field]['datatype']) #,fillvalue=fillvalue)
#            setattr(getattr(self, group), field, this_field[:])
                
 #            for key,val in FH[pt].items():
    
        # if the field dict is not specified, read it
#        if field_dict is None:
#            field_dict={}
#            field_dict['ROOT']=[]
#            for key,val in FH[pt].items():
#                if isinstance(val, h5py.Group):
#                    field_dict[key]=[]
#                    for field in FH[pt][key].keys():
#                        field_dict[key].append(field)
#                if isinstance(val, h5py.Dataset):
#                    field_dict['ROOT'].append(key)
        

#   
    return fileout
    
if __name__=='__main__':
    import argparse
#    parser=argparse.ArgumentParser()
#    parser.add_argument('ATL11_file', type=str)
#    parser.add_argument('--Hemisphere','-H', type=int, default=1, help='1 for Norhtern, -1 for Southern')
#    parser.add_argument('--mosaic', '-m', type=str)
#    parser.add_argument('--out_path', '-o', type=str, help='default is ATL11_file path')
#    parser.add_argument('--pdf', action='store_true', default=False, help='write images to .pdf file')
#    parser.add_argument('--nolog', action='store_true', default=False, help='no writing errors to .log file')
#    args=parser.parse_args()
    fileout = ATL15_write()



