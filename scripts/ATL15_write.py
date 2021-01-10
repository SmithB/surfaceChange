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


def ATL15_write(args):
    
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
            dset.dims[ii].label = scale[dim.strip()]
            if dimScale:
                dset.make_scale(field)
            else:
                if dim.strip().startswith('Nt'):
                    dset.dims[ii].attach_scale(file_obj[scale[dim.strip()]])
                else:
                    if '/' in scale[dim.strip()]:
                        scale_loc=scale[dim.strip()].split('/')[-1]
                    else:
                        scale_loc=scale[dim.strip()]
                    
                    try:
                        dset.dims[ii].attach_scale(group_obj[scale_loc])
                    except Exception as E:
                        print(E)

        for attr in attr_names:
             if 'dimensions' not in attr and 'datatype' not in attr:
                 create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
        if field_attrs[field]['datatype'].startswith('int'):
            dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
        elif field_attrs[field]['datatype'].startswith('float'):
            dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
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
            'file' : ['FH','FH_lag1','FH_lag4'], #,'FH_lag8'],
            'vari' : ['','_lag1','_lag4'], #,'_lag8']
           }
    avgs = ['','_10km','_20km','_40km']

    # establish output file
    kk=0
    fileout = args.output
    if os.path.isfile(fileout):
        os.remove(fileout)
    with h5py.File(fileout.encode('ASCII'),'w') as fo:
        # open data attributes file
        with importlib.resources.path('surfaceChange','resources') as pp:
            with open(os.path.join(pp,'ATL15_output_attrs.csv'),'r', encoding='utf-8-sig') as attrfile:
                reader=list(csv.DictReader(attrfile))
  
        attr_names=[x for x in reader[0].keys() if x != 'field' and x != 'group']

        for kk,ave in enumerate(avgs):
            # field_names = [row['field'] for row in reader if row['group'] == 'height_change'+ave]
            # loop over dz*.h5 files for one ave
            for jj in range(len(lags['file'])):
                filein = args.directory+'/dz'+ave+lags['vari'][jj]+'.h5'
                if not os.path.isfile(filein):
                    print('No file:',args.directory+'/'+os.path.basename(filein))
                    continue
                else:
                    print('Reading file:',args.directory+'/'+os.path.basename(filein))
                lags['file'][jj] = h5py.File(filein,'r')  # file object
                dzg=list(lags['file'][jj].keys())[0]      # dzg is group in input file

                if kk==0:  #establish variables in ROOT
                    for fieldroot in ['year']: 
                        field=fieldroot+lags['vari'][jj]
                        data = np.array(lags['file'][jj][dzg][dz_dict[field]])
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
                        if fieldroot == 'year':
                            make_dataset(field,data,field_attrs,fo,fo,scale,dimScale=True)
    
                if jj==0:  # no lag
                    gh = fo.create_group('height_change'+ave)
                    # spatial dimension scales for the gh
                    for field in ['x','y']:
                        data = np.array(lags['file'][jj][dzg][dz_dict[field]])
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
                        make_dataset(field,data,field_attrs,fo,gh,scale,dimScale=True)
                    
                    for fld in ['cell_area','delta_h','delta_h_sigma']:  # fields that can be ave'd but not lagged
                        field = fld+ave
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
                        if fld in ['cell_area']:
                            data = np.array(lags['file'][jj][dzg][dz_dict[fld]])
                            # fill zeros with invalids
                            data[data==0.0] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max                        
                        if fld in ['delta_h_sigma']:
                            data = np.array(lags['file'][jj][dzg]['sigma_'+dzg])
                        if fld in  ['delta_h']:
                            data = np.array(lags['file'][jj][dzg][dzg])

                        make_dataset(field,data,field_attrs,fo,gh,scale,dimScale=False)

#                    # make output dummy variable, if not in input file                
#                    for field in field_names:
#                        print('field ',field)
#                        if not field.startswith('x') and not field.startswith('y') \
#                        and not field.startswith('delta_h') and 'lag' not in field:
#                            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
#                            print('field_attrs',field_attrs)
#                            if dz_dict.get(field)!=None:
#                                data = np.array(lags['file'][jj][dzg][dz_dict[field]+ave])
#                            else:
#                                # place holder data set for now
#                                dimensions = field_attrs[field]['dimensions'].split(',')
#                                data = np.ndarray(shape=tuple([ii+1 for ii in range(len(dimensions))]),dtype=field_attrs[field]['datatype'])
#
#                            make_dataset(field,data,field_attrs,fo,gh,scale,dimScale=False)
                        
                else:  # one of the lags
                    field = 'dhdt'+lags['vari'][jj]+ave
                    data = np.array(lags['file'][jj][dzg][dzg])
                    field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
                    make_dataset(field,data,field_attrs,fo,gh,scale,dimScale=False)
                    
                    field = 'dhdt'+lags['vari'][jj]+'_sigma'+ave
                    data = np.array(lags['file'][jj][dzg]['sigma_'+dzg])
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
    args=parser.parse_args()
    print('args',args)
    fileout = ATL15_write(args)



