#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 16:40:16 2019

@author: ben
"""

#from ATL11_to_ATL15 
import ATL11_to_ATL15
import numpy as np
import pointCollection as pc
from LSsurf.smooth_xyt_fit import smooth_xyt_fit
from PointDatabase import point_data
import matplotlib.pyplot as plt
import os
from ATL11.plotting_scripts.ATL11_multi_plot import ATL11_multi_plot
import shutil

#index_file='/Volumes/ice2//ben/scf/GL_11/U01/GeoIndex.h5'
index_file='/Volumes/ice2//ben/scf/GL_11/GeoIndex_U01.h5'

xy0=[ -260000., -1860000.]
Wxy=4e4
SRS_proj4='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
E_RMS0={'d2z0_dx2':20000./3000/3000, 'd3z_dx2dt':1000./3000/3000, 'd2z_dxdt':3000/3000, 'd2z_dt2':5000}
#t_span=[2019.25, 2019.75]
t_span=[2018.5, 2019.75]
W={'x':Wxy, 'y':Wxy,'t':np.diff(t_span)}
ctr={'x':xy0[0], 'y':xy0[1], 't':np.mean(t_span)}
spacing={'z0':125., 'dz':2000, 'dt':0.25}
error_spacing={'z0':500,'dz':2000, 'dt':0.25}
reference_epoch=0
N_subset=8
#compute_E=False
max_iterations=6
#mask_file=None
mask_file='/Volumes/ice2/ben/ATL14_test/GimpIceMask_100m_edited.tif'
dzdt_lags=[1, 2]
out_file='/home/ben/temp/ATL15_test_E%d_N%d_loose.h5' % (xy0[0]/1000, xy0[1]/1000)

compute_E=True
if True:
    if True:#False:
        data=ATL11_to_ATL15.read_ATL11(xy0, Wxy, index_file, SRS_proj4)
        data=point_data().from_dict({key:data.__dict__[key] for key in data.fields})
    else:
        data=point_data().from_file(out_file)
    
    S=smooth_xyt_fit(data=data, ctr=ctr, W=W, spacing=spacing, E_RMS=E_RMS0,
                 reference_epoch=reference_epoch, N_subset=N_subset, compute_E=False,
                 bias_params=['rgt','cycle'],  max_iterations=max_iterations, 
                 srs_proj4=SRS_proj4, VERBOSE=True, 
                 mask_file=mask_file, mask_scale={0:10, 1:1}, \
                 dzdt_lags=dzdt_lags);
    ATL11_to_ATL15.save_fit_to_file(S, out_file, dzdt_lags=dzdt_lags, reference_epoch=reference_epoch)
    shutil.copy(out_file, out_file+'_backup')
if compute_E:
    shutil.copy(out_file+'_backup', out_file)
    E_data=point_data().from_file(out_file, group='data')
    SE=smooth_xyt_fit(data=E_data, ctr=ctr, W=W, spacing=error_spacing, E_RMS=E_RMS0,
                 reference_epoch=reference_epoch, N_subset=None, compute_E=True,
                 bias_params=['rgt','cycle'],  max_iterations=0, 
                 srs_proj4=SRS_proj4, VERBOSE=True,
                 mask_file=mask_file, mask_scale={0:10, 1:1}, \
                 dzdt_lags=dzdt_lags, calc_error_file=out_file);
    ATL11_to_ATL15.save_errors_to_file(SE, out_file)




plt.figure(1); plt.clf();  
bounds=[S['grids']['dz'].bds[1], S['grids']['dz'].bds[0]]
D_DEM=pc.grid.data().from_geotif('/Volumes/insar7/ben/ArcticDEM/mosaics/Greenland_100m_v3.0.tif', bounds=bounds)
gx, gy=np.gradient(D_DEM.z, D_DEM.x, D_DEM.y)

extent=np.array(D_DEM.extent)/1000
plt.imshow(gx, extent=extent, vmin=-0.1, vmax=0.1, cmap='gray', origin='lower')
dz_mask=np.ones_like(S['grids']['dz'].mask, dtype=np.float64)
dz_mask[S['grids']['dz'].mask==0]=np.NaN
plt.imshow(S['m']['dzdt_lag1'][:,:,4]*dz_mask, 
           extent=extent, 
           origin='lower', cmap='Spectral', vmin=-5, vmax=5, alpha=0.6)
plt.plot(S['data'].x[S['data'].three_sigma_edit]/1000,S['data'].y[S['data'].three_sigma_edit]/1000,'k.', markersize=0.5)
hb=plt.colorbar()
hb.set_label('dz/dt, m/yr')

if compute_E:
    plt.figure(2); plt.clf();
    plt.imshow(gx, extent=extent, vmin=-0.1, vmax=0.1, cmap='gray', origin='lower')
    plt.imshow(SE['E']['dz'][:,:,2]*dz_mask, 
           extent=extent,
           origin='lower',  alpha=0.6, vmax=0.15)
    plt.plot(S['data'].x[S['data'].three_sigma_edit]/1000,S['data'].y[S['data'].three_sigma_edit]/1000,'w.', markersize=0.5)
    hb=plt.colorbar()
    hb.set_label('dz/dt error, m/yr')
gxm, gym=np.gradient(S['m']['z0'], 250)
plt.figure(3); plt.imshow(gxm, cmap='gray', vmin=-0.1, vmax=0.1, origin='lower', extent=extent)

# code to ID problem tracks:
if False:
    xx=plt.ginput()[0]
    xx=np.round(np.array(xx)/1.e4)*1.e4
    #data=S['data'].copy()
    #xxd=data.x+1j*data.y
    #ii=np.argmin(np.abs(xx[0]+1j*xx[1]-xxd))
    #print(data.rgt[ii])
    files=pc.geoIndex().from_file(index_file).query_xy(xx, get_data=False).keys()
    for file in files:
        plt.figure(); 
        ATL11_multi_plot(os.path.dirname(index_file)+'/'+file)


#S=ATL11_to_ATL15(xy0, Wxy=4e4, ATL11_index=index_file, E_RMS={}, \
#            t_span=[2019.25, 2019.75])