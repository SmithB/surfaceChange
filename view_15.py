#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 20:07:32 2019

@author: ben
"""

import matplotlib.pyplot as plt
import pointCollection as pc
import sys
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
D=[]
for file in sys.argv[1:]:
    D.append({\
              'dz':pc.grid.data().from_h5(file, group='dz/', field_mapping={'z':'dz'}),
              'count':pc.grid.data().from_h5(file, group='dz/', field_mapping={'z':'count'}),
              'z0':pc.grid.data().from_h5(file, group='z0/', field_mapping={'z':'z0'}),
                'data':pc.data().from_h5(file, group='data/')})

plt.clf()
hax=[]
for ii in range(1, 5):
    hax.append(plt.gcf().add_subplot(2, 2, ii))
    
t_slice=[0, 2]

for Di in D:
    dt=Di['dz'].t[t_slice[1]]-Di['dz'].t[t_slice[0]]
    h_im=hax[0].imshow((Di['dz'].z[:,:,t_slice[1]]-Di['dz'].z[:,:,t_slice[0]])/dt, extent=Di['dz'].extent, cmap='Spectral', origin='lower', vmin=-3, vmax=3)
    plt.axis('scaled')
    #temp=Di['data'][Di['data'].three_sigma_edit==1]
    #hax[0].plot(temp.x, temp.y,'.')
    dz0x, dz0y=np.gradient(Di['z0'].z, Di['z0'].x, Di['z0'].y)
    plt.axis('scaled')

    hax[1].imshow(dz0x+dz0y, extent=Di['z0'].extent, cmap='gray', origin='lower', vmin=-0.1, vmax=0.1)
    plt.axis('scaled')

    h_count=hax[2].imshow(np.sum(Di['count'].z>0, axis=2), extent=Di['dz'].extent, origin='lower', vmin=0, vmax=Di['count'].z.shape[2])
    plt.axis('scaled')
    
    ddz=np.std(np.diff(Di['dz'].z, axis=2), axis=2)
    h_ddz=hax[3].imshow(ddz, extent=Di['dz'].extent, origin='lower', vmin=0, vmax=3)
    plt.axis('scaled')

xx=np.array(np.c_[[Di['dz'].x for Di in D]])
XR=[xx.ravel().min(), xx.ravel().max()]
yy=np.array(np.c_[[Di['dz'].y for Di in D]])
YR=[yy.ravel().min(), yy.ravel().max()]


plt.colorbar(mappable=h_im, ax=hax[0])
plt.colorbar(mappable=h_count, ax=hax[2])
for ax in hax:
    ax.set_aspect('equal','datalim')
    ax.set_xlim(XR)
    ax.set_ylim(YR)


plt.show()

def dzplot(xy, D):
    if isinstance(xy,list):
        xy=xy[0]
    x=np.array(xy[0])
    y=np.array(xy[1])
    x.shape=[1,1]
    y.shape=[1,1]
    dzi=[]
    plt.figure()
    row_ind=np.arange(D[0]['dz'].z.shape[0], dtype=int)
    col_ind=np.arange(D[0]['dz'].z.shape[1], dtype=int)

    for Di in D:
        dzsub=[Di['dz'].copy().subset(row_ind, col_ind, band_ind=band).interp(x, y) for band in range(Di['dz'].z.shape[2])]
        dzi.append(dzsub)
        plt.plot(Di['dz'].t, np.array(dzsub).ravel())
        