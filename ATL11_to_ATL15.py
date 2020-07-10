#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 17:15:08 2019

@author: ben
 """
import numpy as np
from LSsurf.smooth_xytb_fit import smooth_xytb_fit
import pointCollection as pc
import os
import re
import sys
import h5py
import matplotlib.pyplot as plt
from surfaceChange.reread_data_from_fits import reread_data_from_fits

def get_SRS_proj4(hemisphere):
    if hemisphere==1:
        return '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
    else:
        return '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

def read_ATL11(xy0, Wxy, index_file, SRS_proj4):
    field_dict_11={'corrected_h':['latitude','longitude','delta_time',\
                        'h_corr','h_corr_sigma','h_corr_sigma_systematic','quality_summary'],\
                        '__calc_internal__':['rgt'],
                        'ref_surf':['e_slope','n_slope', 'x_atc']}
    D11_list=pc.geoIndex().from_file(index_file).query_xy_box(xy0[0]+\
                        np.array([-Wxy/2, Wxy/2]), xy0[1]+np.array([-Wxy/2, Wxy/2]), fields=field_dict_11)
    D_list=[]
    for D11 in D11_list:
        D11.get_xy(proj4_string=SRS_proj4)
        sigma_corr=np.sqrt((7.5*np.abs(np.median(D11.n_slope)))**2+\
                           (7.5*np.abs(np.median(D11.e_slope)))**2+0.03**2)
        # fix for early ATL11 versions that had some zero error estimates.
        bad=np.any(D11.h_corr_sigma==0, axis=1)
        D11.h_corr_sigma[bad,:]=np.NaN

        n_cycles=np.sum(np.isfinite(D11.h_corr), axis=1)
        n_cycles=np.reshape(n_cycles, (D11.shape[0],1))
        n_cycles=np.tile(n_cycles, [1, D11.shape[1]])

        D_list += [pc.data().from_dict({'z':D11.h_corr,
           'sigma_corr':sigma_corr+np.zeros_like(D11.h_corr),
           'sigma':D11.h_corr_sigma,
           'x':D11.x,
           'y':D11.y,
           'x_atc': D11.x_atc,
           'latitude':D11.latitude,
           'longitude':D11.longitude,
           'rgt':D11.rgt,
           'cycle':D11.cycle_number,
           'n_cycles': n_cycles,
           'fit_quality': D11.quality_summary,
           'time':D11.delta_time/24/3600/365.25+2018})]

    D=pc.data().from_list(D_list).ravel_fields()
    D.index(D.fit_quality != 6)
    #bad=(D.rgt==643) & (D.cycle==4)
    #D.index(bad==0)
    return D

def ATL11_to_ATL15(xy0, Wxy=4e4, ATL11_index=None, E_RMS={}, \
            t_span=[2019.25, 2020.5], \
            spacing={'z0':2.5e2, 'dz':5.e2, 'dt':0.25},  \
            dzdt_lags=[1, 4],\
            hemisphere=1, reference_epoch=None, reread_dirs=None, \
            max_iterations=5, \
            N_subset=8,  \
            W_edit=None, \
            out_name=None, replace=False, DOPLOT=False, \
            compute_E=False, mask_file=None, \
            calc_error_file=None):

    SRS_proj4=get_SRS_proj4(hemisphere)
 
    E_RMS0={'d2z0_dx2':200000./3000/3000, 'd3z_dx2dt':3000./3000/3000, 'd2z_dxdt':3000/3000, 'd2z_dt2':5000}
    E_RMS0.update(E_RMS)

    W={'x':Wxy, 'y':Wxy,'t':np.diff(t_span)}
    ctr={'x':xy0[0], 'y':xy0[1], 't':np.mean(t_span)}

    if calc_error_file is None:
        if reread_dirs is None:
            data=read_ATL11(xy0, Wxy, ATL11_index, SRS_proj4)
        else:
            data = reread_data_from_fits(xy0, Wxy, reread_dirs, template='E%d_N%d.h5')
    else:
        data=pc.data().from_h5(calc_error_file, group='data')
        max_iterations=0
        compute_E=True
        N_subset=None
        
    if W_edit is not None:
        # this is used for data that we are rereading from a set of other files.
        # data that are far from the center of this file cannot be eliminated by
        # the three-sigma edit
        data.assign({'editable':  (np.abs(data.x-xy0[0])<=W_edit/2) & (np.abs(data.y-xy0[1])<=W_edit/2)})

    S=smooth_xytb_fit(data=data, ctr=ctr, W=W, spacing=spacing, E_RMS=E_RMS0,
                     reference_epoch=reference_epoch, N_subset=N_subset, compute_E=compute_E,
                     bias_params=['rgt','cycle'],  max_iterations=max_iterations,
                     srs_proj4=SRS_proj4, VERBOSE=True, dzdt_lags=dzdt_lags,
                     mask_file=mask_file, mask_scale={0:10, 1:1})

    #S=smooth_xytb_fit(data=data, ctr=ctr, W=W, spacing=spacing, E_RMS=E_RMS0,
    #                 reference_epoch=reference_epoch, N_subset=N_subset, compute_E=compute_E,
    #                 bias_params=['day','sensor'], repeat_res=repeat_res, max_iterations=max_iterations,
    #                 srs_proj4=SRS_proj4, VERBOSE=True, Edit_only=Edit_only, data_slope_sensors=DEM_sensors,
    #                 mask_file=mask_file, mask_scale={0:10, 1:1})
    return S

def save_fit_to_file(S,  filename, dzdt_lags=None, reference_epoch=0):
    if os.path.isfile(filename):
        os.remove(filename)
    with h5py.File(filename,'w') as h5f:
        h5f.create_group('/data')
        for key in S['data'].fields:
            h5f.create_dataset('/data/'+key, data=getattr(S['data'], key))
        h5f.create_group('/meta')
        h5f.create_group('/meta/timing')
        for key in S['timing']:
            h5f['/meta/timing/'].attrs[key]=S['timing'][key]
            
        h5f.create_group('/dz')
        dz_mask=np.ones_like(S['grids']['dz'].mask, dtype=np.float64)
        dz_mask[S['grids']['dz'].mask==0]=np.NaN
        for ii, name in enumerate(['y','x','t']):
            h5f.create_dataset('/dz/'+name, data=S['grids']['dz'].ctrs[ii])
        sz=list(dz_mask.shape) + [1]
        h5f.create_dataset('/dz/dz', data=S['m']['dz'].dz*np.tile(dz_mask.reshape(sz), [1, 1, S['m']['dz'].shape[2]]))
        h5f.create_dataset('/dz/count', data=S['m']['dz'].count)
        
        h5f.create_group('/z0')
        h5f['z0'].attrs['time']=S['grids']['dz'].ctrs[2][reference_epoch]
        for ii, name in enumerate(['y','x']):
            h5f.create_dataset('/z0/'+name, data=S['grids']['z0'].ctrs[ii])
        z0_mask=np.ones_like(S['grids']['z0'].mask, dtype=np.float64)
        z0_mask[S['grids']['z0'].mask==0]=np.NaN
        h5f.create_dataset('/z0/z0', data=S['m']['z0'].z0*z0_mask)
        
        h5f.create_group('/dz/center_average')
        h5f.create_dataset('/dz/center_average/dz', data=S['m']['dz_bar'])
        for lag in dzdt_lags:
            this_name='dzdt_bar_lag%d' % lag
            h5f.create_dataset('/dz/center_average/'+this_name, data=S['m'][this_name])
        h5f.create_group('/RMS')
        for key in S['RMS']:
            h5f.create_dataset('/RMS/'+key, data=S['RMS'][key])
        h5f.create_group('E_RMS')
        for key in S['E_RMS']:
            h5f.create_dataset('E_RMS/'+key, data=S['E_RMS'][key])
        for key in S['m']['bias']:
            h5f.create_dataset('/bias/'+key, data=S['m']['bias'][key])
    return

def save_errors_to_file( S, filename, dzdt_lags=None, reference_epoch=None):
    
    with h5py.File(filename,'r+') as h5f:
        h5f.create_group('/dz/sigma/')
        for ii, name in enumerate(['y','x','t']):
            h5f.create_dataset('/dz/sigma/'+name, data=S['grids']['dz'].ctrs[ii])
        h5f.create_dataset('dz/sigma/dz', data=S['E']['dz'].dz)
        for lag in dzdt_lags:
            field='dzdt_lag%d' % lag
            h5f.create_dataset('/dz/sigma/'+field, data=S['E'][field])
        h5f.create_group('/z0/sigma')
        for ii, name in enumerate(['y','x']):
            h5f.create_dataset('/z0/sigma/'+name, data=S['grids']['z0'].ctrs[ii])
        h5f.create_dataset('/z0/sigma/z0', data=S['E']['z0'].z0)
        for key in S['E']['bias']:
            h5f.create_dataset('/bias/sigma/'+key, data=S['E']['bias'][key])
       
    return

def main(argv):
    # account for a bug in argparse that misinterprets negative agruents
    for i, arg in enumerate(argv):
        if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

    import argparse
    parser=argparse.ArgumentParser(description="function to fit icebridge data with a smooth elevation-change model", \
                                   fromfile_prefix_chars="@")
    parser.add_argument('xy0', type=float, nargs=2, help="fit center location")
    parser.add_argument('--ATL11_index', type=str, required=True, help="ATL11 index file")
    parser.add_argument('--Width','-W',  type=float, help="Width of grid")
    parser.add_argument('--time_span','-t', type=str, help="time span, first year,last year AD (comma separated, no spaces)")
    parser.add_argument('--grid_spacing','-g', type=str, help='grid spacing:DEM (meters),dh maps xy (meters),dh_maps time (years): comma-separated, no spaces', default='250.,4000.,1.')
    parser.add_argument('--Hemisphere','-H', type=int, default=1, help='hemisphere: -1=Antarctica, 1=Greenland')
    parser.add_argument('--base_directory','-b', type=str, help='base directory')
    parser.add_argument('--out_name', '-o', type=str, help="output file name")
    parser.add_argument('--centers', action="store_true")
    parser.add_argument('--edges', action="store_true")
    parser.add_argument('--corners', action="store_true")
    parser.add_argument('--E_d2zdt2', type=float, default=5000)
    parser.add_argument('--E_d2z0dx2', type=float, default=0.02)
    parser.add_argument('--E_d3zdx2dt', type=float, default=0.0003)
    parser.add_argument('--data_gap_scale', type=float,  default=2500)
    parser.add_argument('--dzdt_lags', type=str, default='1,4', help='lags for which to calculate dz/dt, comma-separated list, no spaces')
    parser.add_argument('--N_subset', type=int, default=None, help="number of pieces into which to divide the domain for (cheap) editing iterations.")
    parser.add_argument('--max_iterations', type=int, default=6, help="maximum number of iterations used to edit the data.")
    parser.add_argument('--map_dir','-m', type=str)
    parser.add_argument('--mask_file', type=str)
    parser.add_argument('--reference_epoch', type=int, default=0, help="Reference epoch number, for which dz=0")
    parser.add_argument('--calc_error_file','-c', type=str, help='file containing data for which errors will be calculated')
    parser.add_argument('--calc_error_for_xy', action='store_true', help='calculate the errors for the file specified by the x0, y0 arguments')
    parser.add_argument('--error_res_scale','-s', type=float, default=2, help='if the errors are being calculated (see calc_error_file), scale the grid resolution in x and y to be coarser')
    args=parser.parse_args()

    args.grid_spacing = [np.float(temp) for temp in args.grid_spacing.split(',')]
    args.dzdt_lags = [np.int(temp) for temp in args.dzdt_lags.split(',')]
    args.time_span = [np.float(temp) for temp in args.time_span.split(',')]

    spacing={'z0':args.grid_spacing[0], 'dz':args.grid_spacing[1], 'dt':args.grid_spacing[2]}
    E_RMS={'d2z0_dx2':args.E_d2z0dx2, 'd3z_dx2dt':args.E_d3zdx2dt, 'd2z_dxdt':args.E_d3zdx2dt*args.data_gap_scale,  'd2z_dt2':args.E_d2zdt2}
    
    reread_dirs=None
    dest_dir=args.base_directory
    W_edit=args.Width/2
    if args.centers:
        dest_dir += '/centers'
        W_edit=None
    if args.edges or args.corners:
        reread_dirs=[args.base_directory+'/centers']
        if args.edges:
            dest_dir += '/edges'
    if args.corners:
        reread_dirs += [args.base_directory+'/edges']
        dest_dir +='/corners'

    if args.out_name is None:
        args.out_name=dest_dir + '/E%d_N%d.h5' % (args.xy0[0]/1e3, args.xy0[1]/1e3)

    if args.calc_error_file is not None:
        dest_dir=os.path.dirname(args.calc_error_file)
        # get xy0 from the filename
        re_match=re.compile('E(.*)_N(.*).h5').search(args.calc_error_file)
        args.xy0=[float(re_match.group(ii))*1000 for ii in [1, 2]]

    if args.calc_errors_for_xy is not None:
        args.calc_error_file=args.out_name

    if args.error_res_scale is not None:
        if args.calc_error_file is None:
            print("if error_res_scale is specified, either --calc_errors_for_xy  or --calc_error_file must be present")
            sys.exit(1)

    if not os.path.isdir(args.base_directory):
        os.mkdir(args.base_directory)
    if not os.path.isdir(dest_dir):
        os.mkdir(dest_dir)

    S=ATL11_to_ATL15(args.xy0, ATL11_index=args.ATL11_index,
           Wxy=args.Width, E_RMS=E_RMS, t_span=args.time_span, spacing=spacing, \
           hemisphere=args.Hemisphere, reread_dirs=reread_dirs, \
           out_name=args.out_name, 
           DOPLOT=False, \
           dzdt_lags=args.dzdt_lags, \
           N_subset=args.N_subset,\
           mask_file=args.mask_file, \
           max_iterations=args.max_iterations, \
           reference_epoch=args.reference_epoch, \
           W_edit=W_edit,
           calc_error_file=args.calc_error_file)
    
    if args.calc_error_file is None:
        # if this isn't an error-calculation run, save the gridded fit data to the output file
        save_fit_to_file(S, args.out_name, dzdt_lags=args.dzdt_lags, reference_epoch=args.reference_epoch)
    else:
        # If this is an error-calculation run, save the errors to the output file
        save_errors_to_file(S, args.out_name, dzdt_lags=args.dzdt_lags, reference_epoch=args.reference_epoch)
    print(f"done with {args.out_name}")
if __name__=='__main__':
    main(sys.argv)

#-160000 -1800000 --centers @/home/ben/git_repos/surfaceChange/default_args_z03xlooser_dt10xlooser.txt
#-160000 -1800000 --centers @/home/ben/git_repos/surfaceChange/default_args_z03xlooser_dt10xlooser_errors.txt -c /Volumes/ice2/ben/ATL14_test/IS2//U07/z03xlooser_dt10xlooser_40km/centers/E-160_N-1800.h5