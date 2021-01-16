#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 08:37:27 2020

@author: ben
"""


from LSsurf.fd_grid import fd_grid
from LSsurf.lin_op import lin_op
import scipy.sparse as sp
import matplotlib.pyplot as plt
import numpy as np
from LSsurf.smooth_xytb_fit import sum_cell_area
from LSsurf.smooth_xytb_fit import calc_cell_area
from LSsurf.smooth_xytb_fit import setup_averaging_ops


xc=np.array([0, -5.e5])
deltas=[100., 100.]
bounds=[xc[1]+np.array([-3.e4, 3.e4]), xc[0]+np.array([-3.e4, 3.e4])]
srs_proj4='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '

grid_z0=fd_grid(bounds, deltas, srs_proj4=srs_proj4)
grid_dz=fd_grid(bounds+[np.array([0, 5])], [1.e3, 1.e3, 0.25])
grid_10km=fd_grid(bounds+[np.array([0, 5])], [1.e4, 1.e4, 0.25])

mask=np.zeros(grid_z0.ctrs[0].size*np.array([1,1]))

mask[:,grid_z0.ctrs[1]<np.mean(grid_z0.ctrs[1])-150]=1
grid_z0.mask=mask

cell_area_0 = calc_cell_area(grid_z0)
cell_area_1, op = sum_cell_area(grid_z0, grid_dz, return_op=True)

args={'avg_scales':[1.e4], 'dzdt_lags':[1, 4]}

ops=setup_averaging_ops(grid_dz, grid_dz.col_N, args, cell_area=cell_area_1)


if False:
    temp=op.toCSR().dot(mask.ravel())
    fig=plt.figure();
    fig.add_subplot(131); plt.imshow(cell_area_0*mask, origin='lower')
    fig.add_subplot(132); plt.imshow(op.toCSR()[6,:].toarray().reshape(mask.shape), origin='lower')
    fig.add_subplot(133); plt.imshow(cell_area_1, origin='lower')
