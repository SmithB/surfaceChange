#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 17:15:08 2019

@author: ben
 """
import os
#os.environ["MKL_NUM_THREADS"]="1"  # multiple threads don't help that much and tend to eat resources
os.environ["OPENBLAS_NUM_THREADS"]="1"  # multiple threads don't help that much and tend to eat resources

import numpy as np
from LSsurf.smooth_xytb_fit import smooth_xytb_fit
import pointCollection as pc

import re
import sys
import h5py
import matplotlib.pyplot as plt
from surfaceChange.reread_data_from_fits import reread_data_from_fits
from pyTMD import compute_tide_corrections

import pdb

def get_SRS_proj4(hemisphere):
    if hemisphere==1:
        return '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
    else:
        return '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

def manual_edits(D):
    '''
    Remove known problematic data from a data structure
    
    inputs:
        D: data structure
    outputs:
        None (modifies input structure in place)
    '''
    bad=(D.rgt == 741)  & (D.cycle == 7)
    D=D[~bad]
    return

def read_ATL11(xy0, Wxy, index_file, SRS_proj4):
    '''
    read ATL11 data from an index file
    
    inputs:
        xy0 : 2-element iterable specifying the domain center
        Wxy : Width of the domain
        index_file : file made by pointCollection.geoindex pointing at ATL11 data
        SRS_proj4: projection information for the data

    output:
        D: data structure
    '''
    
    field_dict_11={None:['latitude','longitude','delta_time',\
                        'h_corr','h_corr_sigma','h_corr_sigma_systematic', 'ref_pt'],\
                        '__calc_internal__' : ['rgt'],
                        'cycle_stats' : {'tide_ocean','dac'},
                        'ref_surf':['e_slope','n_slope', 'x_atc', 'fit_quality']}
    try:
        # catch empty data
        D11_list=pc.geoIndex().from_file(index_file).query_xy_box(
            xy0[0]+np.array([-Wxy/2, Wxy/2]), \
            xy0[1]+np.array([-Wxy/2, Wxy/2]), fields=field_dict_11)
    except ValueError:
        return None
    D_list=[]
    XO_list=[]
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
           'ref_pt':D11.ref_pt,
           'cycle':D11.cycle_number,
           'n_cycles': n_cycles,
           'fit_quality': D11.fit_quality,
           'tide_ocean': D11.tide_ocean,
           'dac': D11.dac,
           'delta_time': D11.delta_time,
           'time':D11.delta_time/24/3600/365.25+2018})]
        if len(D11.ref_pt) == 0:
            continue
        # N.B.  D11 is getting indexed in this step, and it's leading to the warning in 
        # line 76.  Can fix by making crossover_data.from_h5 copy D11 on input
        D_x = pc.ATL11.crossover_data().from_h5(D11.filename, pair=D11.pair, D_at=D11)
        if D_x is None:
            continue
        # fit_quality D11 applies to all cycles, but is mapped only to some specific
        # cycle.
        temp=np.nanmax(D_x.fit_quality[:,:,0], axis=1)
        for cycle in range(D_x.shape[1]): 
            D_x.fit_quality[:,cycle,1]=temp

        D_x.get_xy(proj4_string=SRS_proj4)

        good=np.isfinite(D_x.h_corr)[:,0:2,1].ravel()
        for field in D_x.fields:
            # Pull out only cycles 1 and 2
            temp=getattr(D_x, field)[:,0:2,1]
            setattr(D_x, field, temp.ravel()[good])

        zero = np.zeros_like(D_x.h_corr)
        blank = zero+np.NaN
        XO_list += [pc.data().from_dict({'z':D_x.h_corr,
            'sigma':D_x.h_corr_sigma,
            'sigma_corr': sigma_corr+np.zeros_like(D_x.h_corr),
            'x':D_x.x,
            'y':D_x.y,
            'latitude':D_x.latitude,
            'longitude':D_x.longitude,
            'rgt':D_x.rgt,
            'ref_pt':blank,
            'cycle':D_x.cycle_number,
            'n_cycles':blank,
            'fit_quality':D_x.fit_quality,
            'tide_ocean':D_x.tide_ocean,
            'dac':D_x.dac,
            'delta_time':D_x.delta_time,
            'time':D_x.delta_time/24/3600/365.25+2018})]
    try:
        D=pc.data().from_list(D_list+XO_list).ravel_fields()
    except ValueError:
        # catch empty data
        return None
    D.index( ( D.fit_quality ==0 ) | ( D.fit_quality == 2 ))
    
    return D

def apply_tides(D, xy0, W, tide_mask_file, tide_directory):
    '''
    read in the tide mask (for Antarctica) and apply dac and tide to ice-shelf elements
    '''
    # the tide mask should be 1 for non-grounded points (ice shelves?), zero otherwise
    tide_mask = pc.grid.data().from_geotif(tide_mask_file, bounds=[np.array([-0.6, 0.6])*W+xy0[0], np.array([-0.6, 0.6])*W+xy0[1]])     
    is_els=tide_mask.interp(D.x, D.y) > 0.5
    print(f"\t\t{np.mean(is_els)*100} shelf data")
    if np.any(is_els.ravel()):
        D.tide_ocean = compute_tide_corrections(\
                D.x, D.y, D.delta_time,                                                           
                DIRECTORY=tide_directory, MODEL='CATS2008',                                            
                EPOCH=(2018,1,1,0,0,0), TYPE='drift', TIME='utc', EPSG=3031)    
    D.tide_ocean[is_els==0] = 0
    D.dac[is_els==0] = 0
    D.tide_ocean[~np.isfinite(D.tide_ocean)] = 0
    D.dac[~np.isfinite(D.dac)] = 0
    D.z -= (D.tide_ocean + D.dac)
    return D

def decimate_data(D, N_target, W_domain,  W_sub, x0, y0):
    # reduce the data density to a target value in small bins
    ij_bin=np.round((D.x - (x0 - W_domain/2 + W_sub/2))/W_sub)+ \
        1j*np.round((D.y - (y0 - W_domain/2 + W_sub/2))/W_sub)
    rho_target = N_target / W_domain**2
    ind_buffer=[]
    for bin0 in np.unique(ij_bin):
        ii = np.flatnonzero(ij_bin==bin0)
        this_rho = len(ii) / (W_sub**2)
        #print(f'bin0={bin0}, this_rho={this_rho}, ratio={this_rho/rho_target}')
        if this_rho < rho_target:
            # if there aren't too many points in this bin, continue
            ind_buffer += [ii] 
            continue
        # make a global reference point number (equal to the number of ref pts in an orbit * rgt + ref_pt)
        global_ref_pt=D.ref_pt[ii]+40130000/20*D.rgt[ii]
        u_ref_pts = np.unique(global_ref_pt)
        # select a subset of the global reference points (this skips the right number of points)
        # note that the np.arange() call can sometimes return a value of len(u_ref_pts) (?), so 
        # its output needs to be checked (!)
        sel_ind=np.arange(0, len(u_ref_pts), this_rho/rho_target).astype(int)
        sel_ref_pts = u_ref_pts[sel_ind[sel_ind < len(u_ref_pts)]]       
        isub = np.in1d(global_ref_pt, sel_ref_pts)
        ind_buffer.append(ii[isub])
    D.index(np.concatenate(ind_buffer))
    
def ATL11_to_ATL15(xy0, Wxy=4e4, ATL11_index=None, E_RMS={}, \
            t_span=[2019.25, 2020.5], \
            spacing={'z0':2.5e2, 'dz':5.e2, 'dt':0.25},  \
            dzdt_lags=[1, 4],\
            hemisphere=1, reference_epoch=None, reread_dirs=None, \
            data_file=None, \
            max_iterations=5, \
            N_subset=8,  \
            W_edit=None, \
            out_name=None, \
            compute_E=False, \
            mask_file=None, \
            tide_mask_file=None,\
            tide_directory=None,\
            avg_scales=None,\
            edge_pad=None,\
            error_res_scale=None,\
            calc_error_file=None):
    '''
    Function to generate DEMs and height-change maps based on ATL11 surface height data.
    
    Inputs:
        xy0: 2 elements specifying the center of the domain
        Wxy: (float) Width of the domain
        ATL11_index: (string) Index file (from pointCollection.geoIndex) for ATL11 data
        E_RMS: (dict) dictionary specifying constraints on derivatives of the ice-sheet surface.  
        t_span: (2-element list of floats) starting and ending times for the output grids (in years CE)
        spacing: (dict) dictionary specifying the grid spacing for z0, dz, and dt
        dzdt_lags: (list) lags over which elevation change rates and errors will be calculated
        hemisphere: (int) the hemisphere in which the grids will be generated. 1 for northern hemisphere, -1 for southern
        reference_epoch: (int) The time slice (counting from zero, in steps of spacing['dt']) corresponding to the DEM
        reread_dirs: (string) directory containing output files from which data will be read (if None, data will be read from the index file)
        data_file: (string) output file from which to reread data (alternative to reread_dirs)
        max_iterations: (int) maximum number of iterations in three-sigma edit of the solution
        N_subset: (int) If specified, the domain is subdivided into this number of divisions in x and y, and a three-sigma edit is calculated for each
        W_edit: (float) Only data within this distance of the grid center can be editied in the iterative step
        out_name: (string) Name of the output file
        compute_E: (bool) If true, errors will be calculated
        mask_file: (string) File specifying areas for which data should be used and strong constraints should be applied
        tide_mask_file: (string)  File specifying the areas for which the tide model should be calculated
        tide_directory: (string)  Directory containing the tide model data
        avg_scales: (list of floats) scales over which the output grids will be averaged and errors will be calculated
        error_res_scale: (float) If errors are calculated, the grid resolution will be coarsened by this factor
        calc_error_file: (string) Output file for which errors will be calculated.
    '''
    SRS_proj4=get_SRS_proj4(hemisphere)
 
    E_RMS0={'d2z0_dx2':200000./3000/3000, 'd3z_dx2dt':3000./3000/3000, 'd2z_dxdt':3000/3000, 'd2z_dt2':5000}
    E_RMS0.update(E_RMS)

    W={'x':Wxy, 'y':Wxy,'t':np.diff(t_span)}
    ctr={'x':xy0[0], 'y':xy0[1], 't':np.mean(t_span)}

    # figure out where to get the data
    if data_file is not None:
        data=pc.data().from_h5(data_file, group='/')
    elif calc_error_file is not None:
        data=pc.data().from_h5(calc_error_file, group='data')
        max_iterations=0
        compute_E=True
        N_subset=None
    elif reread_dirs is not None:
        data = reread_data_from_fits(xy0, Wxy, reread_dirs, template='E%d_N%d.h5')
    else:
        data=read_ATL11(xy0, Wxy, ATL11_index, SRS_proj4)
        N0=data.size
        decimate_data(data, 1.2e6, Wxy, 5000, xy0[0], xy0[1] )
        print(f'decimated {N0} to {data.size}')
    
    if data is None:
        print("No data present for region, returning.")
        return None
    
    # if any manual edits are needed, make them here:
    manual_edits(data)
    
    if data.time.max()-data.time.min() < 80./365.:
        print("time span too short, returning.")
        return None
    
    if edge_pad is not None:
        ctr_dist = np.max(np.abs(data.x-xy0[0]), np.abs(data.y-xy0[1]))
        in_ctr = ctr_dist < Wxy/2 - edge_pad
        if np.sum(in_ctr) < 50:
            return None
    
    # apply the tides if a directory has been provided
    if tide_mask_file is not None:
        apply_tides(data, xy0, Wxy, tide_mask_file, tide_directory)
        
    if W_edit is not None:
        # this is used for data that we are rereading from a set of other files.
        # data that are far from the center of this file cannot be eliminated by
        # the three-sigma edit
        data.assign({'editable':  (np.abs(data.x-xy0[0])<=W_edit/2) & (np.abs(data.y-xy0[1])<=W_edit/2)})
    
    # call smooth_xytb_fitting
    S=smooth_xytb_fit(data=data, ctr=ctr, W=W, spacing=spacing, E_RMS=E_RMS0,
                     reference_epoch=reference_epoch, N_subset=N_subset, compute_E=compute_E,
                     bias_params=['rgt','cycle'],  max_iterations=max_iterations,
                     srs_proj4=SRS_proj4, VERBOSE=True, dzdt_lags=dzdt_lags,
                     mask_file=mask_file, mask_scale={0:10, 1:1}, 
                     error_res_scale=error_res_scale,
                     avg_scales=avg_scales)
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
        h5f.create_group('/RMS')
        for key in S['RMS']:
           h5f.create_dataset('/RMS/'+key, data=S['RMS'][key])
        h5f.create_group('E_RMS')
        for key in S['E_RMS']:
            h5f.create_dataset('E_RMS/'+key, data=S['E_RMS'][key])
        for key in S['m']['bias']:
            h5f.create_dataset('/bias/'+key, data=S['m']['bias'][key])
    for key , ds in S['m'].items():
        if isinstance(ds, pc.grid.data):
            ds.to_h5(filename, group=key)
    return

def interp_ds(ds, scale):
    for field in ds.fields:
        delta_xy=[(ds.x[1]-ds.x[0])/scale, (ds.y[1]-ds.y[0])/scale]
        xi=np.arange(ds.x[0], ds.x[-1]+delta_xy[0], delta_xy[0])
        yi=np.arange(ds.y[0], ds.y[-1]+delta_xy[1], delta_xy[1])
        z0=getattr(ds, field)
        if len(ds.shape)==2:
            zi=pc.grid.data().from_dict({'x':ds.x, 'y':ds.y, 'z':z0}).interp(xi, yi, gridded=True)
            return pc.grid.data().from_dict({'x':xi, 'y':yi, field:zi})
        else:
            zi=np.zeros([xi.size, yi.size, ds.time.size])
            for epoch in range(ds.time.size):
                temp=pc.grid.data().from_dict({'x':ds.x, 'y':ds.y, 'z':np.squeeze(z0[:,:,epoch])})
                zi[:,:,epoch] = temp.interp(xi, yi, gridded=True)
            return pc.grid.data().from_dict({'x':xi, 'y':yi, 'time':ds.time, field:zi})

def save_errors_to_file( S, filename, dzdt_lags=None, reference_epoch=None, grid_datasets=None):

    for key, ds in S['E'].items():
        if isinstance(ds, pc.grid.data):
            print(key)
            ds.to_h5(filename, group=key.replace('sigma_',''))

    with h5py.File(filename,'r+') as h5f:
        for key in S['E']['sigma_bias']:
            if 'bias/sigma' in h5f and  key in h5f['/bias/sigma']:
                print(f'{key} already exists in sigma_bias')
                h5f['/bias/sigma/'+key][...]=S['E']['sigma_bias'][key]
            else:
                h5f.create_dataset('/bias/sigma/'+key, data=S['E']['sigma_bias'][key])
       
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
    parser.add_argument('--avg_scales', type=str, default='10000,40000', help='scales at which to report average errors, comma-separated list, no spaces')
    parser.add_argument('--N_subset', type=int, default=None, help="number of pieces into which to divide the domain for (cheap) editing iterations.")
    parser.add_argument('--max_iterations', type=int, default=6, help="maximum number of iterations used to edit the data.")
    parser.add_argument('--map_dir','-m', type=str)
    parser.add_argument('--mask_file', type=str)
    parser.add_argument('--tide_mask_file', type=str)
    parser.add_argument('--tide_directory', type=str)
    parser.add_argument('--reference_epoch', type=int, default=0, help="Reference epoch number, for which dz=0")
    parser.add_argument('--data_file', type=str, help='read data from this file alone')
    parser.add_argument('--calc_error_file','-c', type=str, help='file containing data for which errors will be calculated')
    parser.add_argument('--calc_error_for_xy', action='store_true', help='calculate the errors for the file specified by the x0, y0 arguments')
    parser.add_argument('--error_res_scale','-s', type=float, nargs=2, default=[4, 2], help='if the errors are being calculated (see calc_error_file), scale the grid resolution in x and y to be coarser')
    args=parser.parse_args()


    args.grid_spacing = [np.float(temp) for temp in args.grid_spacing.split(',')]
    args.dzdt_lags = [np.int(temp) for temp in args.dzdt_lags.split(',')]
    args.time_span = [np.float(temp) for temp in args.time_span.split(',')]
    args.avg_scales = [np.int(temp) for temp in args.avg_scales.split(',')]

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
        args.out_name=args.calc_error_file
        if not os.path.isfile(args.out_name):
            print(f"{args.out_name} not found, returning")
            return
        
    if args.calc_error_for_xy:
        args.calc_error_file=args.out_name
        if not os.path.isfile(args.out_name):
            print(f"{args.out_name} not found, returning")
            return
        
    if args.error_res_scale is not None:
        if args.calc_error_file is not None:
            for ii, key in enumerate(['z0','dz']):
                spacing[key] *= args.error_res_scale[ii]

    if not os.path.isdir(args.base_directory):
        os.mkdir(args.base_directory)
    if not os.path.isdir(dest_dir):
        os.mkdir(dest_dir)

    S=ATL11_to_ATL15(args.xy0, ATL11_index=args.ATL11_index,
           Wxy=args.Width, E_RMS=E_RMS, t_span=args.time_span, spacing=spacing, \
           hemisphere=args.Hemisphere, reread_dirs=reread_dirs, \
           data_file=args.data_file, \
           out_name=args.out_name,
           dzdt_lags=args.dzdt_lags, \
           N_subset=args.N_subset,\
           mask_file=args.mask_file, \
           tide_mask_file=args.tide_mask_file, \
           tide_directory=args.tide_directory, \
           max_iterations=args.max_iterations, \
           reference_epoch=args.reference_epoch, \
           W_edit=W_edit,\
           calc_error_file=args.calc_error_file, \
           error_res_scale=args.error_res_scale, \
           avg_scales=args.avg_scales)

    if S is None:
        return
    
    if args.calc_error_file is None and 'm' in S and len(S['m'].keys()) > 0:
        # if this isn't an error-calculation run, save the gridded fit data to the output file
        save_fit_to_file(S, args.out_name, dzdt_lags=args.dzdt_lags, reference_epoch=args.reference_epoch)
    elif 'E' in S and len(S['E'].keys()) > 0:
        # If this is an error-calculation run, save the errors to the output file
        S['E']['sigma_z0']=interp_ds(S['E']['sigma_z0'], args.error_res_scale[0])
        for field in ['sigma_dz'] + [ f'sigma_dzdt_lag{lag}' for lag in args.dzdt_lags ]:
            S['E'][field] = interp_ds( S['E'][field], args.error_res_scale[1] )
        save_errors_to_file(S, args.out_name, dzdt_lags=args.dzdt_lags, reference_epoch=args.reference_epoch)
    print(f"done with {args.out_name}")
if __name__=='__main__':
    main(sys.argv)

#-160000 -1800000 --centers @/home/ben/git_repos/surfaceChange/default_args/test.txt
#-160000 -1800000 --centers @/home/ben/git_repos/surfaceChange/default_args_z03xlooser_dt10xlooser_errors.txt -c /Volumes/ice2/ben/ATL14_test/IS2//U07/z03xlooser_dt10xlooser_40km/centers/E-160_N-1800.h5
