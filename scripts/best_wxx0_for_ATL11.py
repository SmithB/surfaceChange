#! /bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from LSsurf import smooth_xytb_fit
import pointCollection as pc
#import sparseqr
#import glob
import h5py
#import os
import LSsurf
#import re
import sys

def safe_interp(x, x0_in, y0_in):
    y=np.NaN
    
    if x0_in[-1] < x0_in[0]:
        x0=x0_in[::-1]
        y0=y0_in[::-1]
    else:
        x0=x0_in
        y0=y0_in
    try:
        i0=np.argwhere(x0 < x)[-1][0]
        i1=np.argwhere(x0 >=x)[0][0]
        #print([i0, i1])
        #print( x0[[i0, i1]])
        #print( y0[[i0, i1]])
        y=np.interp(x, x0[[i0, i1]], y0[[i0, i1]])
    except Exception:
        pass
    return y

def read_ATL11_file(file, mask_file):

    D11=pc.ATL11.data().from_h5(file)
    with h5py.File(file,'r') as h5f:
        qs=np.array(h5f['/pt2/ref_surf/fit_quality'])
    D11.assign({'ref_surf_quality':qs})
    D11.get_xy(EPSG=3413)
    XR=np.array([np.nanmin(D11.x), np.nanmax(D11.x)])
    YR=np.array([np.nanmin(D11.y), np.nanmax(D11.y)])
    mask=pc.grid.data().from_geotif(mask_file, bounds=[XR, YR]).interp(D11.x[:,0], D11.y[:,0]) > 0.5
    D11.index(mask & (D11.ref_surf_quality <1))
    return D11

def save_L_interp(L_interp, file):
    pt0=np.array(list(L_interp.keys()))
    D={}
    for field in L_interp[pt0[0]]:
        D[field]=np.array([L_interp[pt][field] for pt in pt0])
    pc.data().from_dict(D).to_h5(file)


def find_best_wxx0(D11, mask_file):

    scale_vals=np.array([ 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300])

    E_d3zdx2dt=0.0001
    E_d2z0dx2=0.006
    E_d2zdt2=5000

    data_gap_scale=2500

    # define the domain's width in x, y, and time
    W={'x':4.e4,'y':200,'t':.2}
    # define the grid center:
    XR=np.nanmean(D11.x_atc)+np.array([-1, 1])*W['x']/2
    ctr={'x':XR[0]+W['x']/2., 'y':0., 't':0.}
    # define the grid spacing
    spacing={'z0':100, 'dz':100, 'dt':.1}

    dN=np.ceil(W['x']/30).astype(int)

    L_interp={}
    for pt0 in np.arange(D11.ref_pt[0,0]+dN/2, D11.ref_pt[-1,0], dN):
        ii=np.flatnonzero(np.abs(D11.ref_pt[:,0]-pt0)<3*dN/2)
        N_good=np.sum(np.isfinite(D11.h_corr[ii,:]), axis=0)
        if np.max(N_good)<0.9*dN:
            continue
        bc=np.argmax(N_good)
        nb=N_good[bc]

        xy_ctr=[np.nanmean(D11.x[ii, bc]), np.nanmean(D11.y[ii, bc]), np.nanmean(D11.h_corr[ii, bc])]

        D=pc.data().from_dict({'x':D11.x_atc[ii,bc], 'y':np.zeros_like(ii, dtype=float),'z':D11.h_corr[ii,bc],\
                           'time':np.zeros_like(ii, dtype=float), 'sigma':D11.h_corr_sigma[ii,bc]})
        D.index(np.isfinite(D.z) & np.isfinite(D.sigma) & (D.sigma>0))
        S=[]
        ctr={'x':np.nanmean(D.x), 'y':0., 't':0.}

        L_curve={key:[] for key in ['wzz0', 'sigma_hat_s', 'N']}

        for  scale_val in scale_vals:
            # run the fit
            E_RMS={'d2z0_dx2': E_d2z0dx2*scale_val,
             'dz0_dx': E_d2z0dx2*data_gap_scale*scale_val,
             'd3z_dx2dt':E_d3zdx2dt ,
             'd2z_dxdt': E_d3zdx2dt*data_gap_scale,
             'd2z_dt2': E_d2zdt2}
            try:
                S.append(smooth_xytb_fit(data=D, ctr=ctr, W=W, spacing=spacing, E_RMS=E_RMS,
                             reference_epoch=1, N_subset=None, compute_E=False,
                             max_iterations=5,
                             VERBOSE=False))
            except Exception as e:
                S.append(smooth_xytb_fit(data=D, ctr=ctr, W=W, spacing=spacing, E_RMS=E_RMS,
                             reference_epoch=1, N_subset=None, compute_E=False,
                             max_iterations=5,
                             VERBOSE=False))
            d_ed = S[-1]['data']
            d_ed.index(d_ed.three_sigma_edit==1)

            L_curve['sigma_hat_s'].append( LSsurf.RDE((d_ed.z-d_ed.z_est)/d_ed.sigma))
            L_curve['wzz0'].append(E_RMS['d2z0_dx2'])
            L_curve['N'].append(d_ed.size)
        for key in L_curve.keys():
            L_curve[key] = np.array(L_curve[key])

        L_interp[pt0] = {"w_for_r_of_1":safe_interp(1, L_curve['sigma_hat_s'], L_curve['wzz0']), 
                        'w_for_r_10pct_above_min':safe_interp(1.1*L_curve['sigma_hat_s'].min(), L_curve['sigma_hat_s'], L_curve['wzz0']), 
                        'x': xy_ctr[0], 
                        'y': xy_ctr[1],
                        'z': xy_ctr[2]}
    return L_interp


ATL11_file=sys.argv[1]
mask_file=sys.argv[2]
out_file=sys.argv[3]
D11=read_ATL11_file(ATL11_file, mask_file)
if D11.size > 0:
    L_interp=find_best_wxx0(D11, mask_file)
if len(list(L_interp.keys())) > 0:
    save_L_interp(L_interp, out_file)
