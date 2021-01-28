#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:45:47 2020

@author: ben05
"""
import numpy as np
import  os, h5py,  csv
import importlib.resources
from ATL11.h5util import create_attribute
from netCDF4 import Dataset
import alphashape


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

        if field.startswith('year'):
            dsetvar = file_obj.createVariable(field,
                                              nctype[field_attrs[field]['datatype']],
                                              dimensions,
                                              fill_value=fill_value)
        else:
            dsetvar = group_obj.createVariable(field,
                                               nctype[field_attrs[field]['datatype']],
                                               dimensions,
                                               fill_value=fill_value)
            
        dsetvar[:] = data
        for attr in attr_names:
            dsetvar.setncattr(attr,field_attrs[field][attr])

        return file_obj
    
    dz_dict ={'year':'t',     # for non-lagged vars. {ATL15 outgoing var name: hdf5 incoming var name}          
              'year_lag1':'t',
              'year_lag4':'t',
              'x':'x',
              'y':'y',
              'cell_area':'cell_area',
              'ice_mask':'mask',
              'data_count':'count',
              'misfit_rms':'misfit_rms',
              'misfit_scaled_rms':'misfit_scaled_rms',
              'delta_h':'dz',
              'delta_h_sigma':'sigma_dz',
              }
    nctype = {'float64':'f8',
              'float32':'f4',
              'int8':'i1'}
    scale = {'Nt':'year',
             'Nt_lag1':'year_lag1',
             'Nt_lag4':'year_lag4',
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

                if kk==0:  #establish variables in ROOT
                    for fieldroot in ['year']: 
                        field=fieldroot+lags['vari'][jj]
                        data = np.array(lags['file'][jj][dzg][dz_dict[field]])
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='ROOT'}
                        if fieldroot == 'year':
                            make_dataset(field,data,field_attrs,nc,nc,scale,nctype,dimScale=True)
    
                if jj==0:  # no lag
                    nc.createGroup('height_change'+ave)
                    # spatial dimension scales for the group
                    for field in ['x','y']:
                        data = np.array(lags['file'][jj][dzg][dz_dict[field]])
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                        make_dataset(field,data,field_attrs,nc,nc.groups['height_change'+ave],scale,nctype,dimScale=True)
                    
                    for fld in ['cell_area','delta_h','delta_h_sigma','misfit_rms','misfit_scaled_rms']:  # fields that can be ave'd but not lagged
                        if kk>0 and fld.startswith('misfit'):
                            break
                        field = fld+ave
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                        if fld in ['cell_area','misfit_rms','misfit_scaled_rms']:
                            data = np.array(lags['file'][jj][dzg][dz_dict[fld]])
                            data = np.moveaxis(data,0,1)
                            # fill zeros with invalids
                            data[data==0.0] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max                        
                        if fld in ['delta_h_sigma']:
                            data = np.array(lags['file'][jj][dzg]['sigma_'+dzg])
                            data = np.moveaxis(data,0,1)
                        if fld in  ['delta_h']:
                            data = np.array(lags['file'][jj][dzg][dzg])
                            data = np.moveaxis(data,0,1)
                        make_dataset(field,data,field_attrs,nc,nc.groups['height_change'+ave],scale,nctype,dimScale=False)

                else:  # one of the lags
                    field = 'dhdt'+lags['vari'][jj]+ave
                    data = np.array(lags['file'][jj][dzg][dzg])
                    data = np.moveaxis(data,0,1)
                    field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                    make_dataset(field,data,field_attrs,nc,nc.groups['height_change'+ave],scale,nctype,dimScale=False)
                    
                    field = 'dhdt'+lags['vari'][jj]+'_sigma'+ave
                    data = np.array(lags['file'][jj][dzg]['sigma_'+dzg])
                    data = np.moveaxis(data,0,1)
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



