#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 12:55:15 2019

@author: ben
"""

import h5py
import matplotlib.pyplot as plt
import pointCollection as pc
import numpy as np
files=['ATL15_test_E-260_N-1860_loose.h5', 'ATL15_test_E-260_N-1860_tight.h5']

DD=[];

for file in files:
    with h5py.File('/home/ben/temp/'+file,'r') as h5f:
        D=dict()
        for field in ['dz','dzdt_lag1','dzdt_lag2']:
            D[field]=pc.grid.data().from_dict(\
             {'x':np.array(h5f['dz']['x']), \
             'y':np.array(h5f['dz']['y']), \
             'time':np.array(h5f['dz']['t']), \
             'z':np.array(h5f['dz'][field])})
            
            D['sigma_'+field]=pc.grid.data().from_dict(\
             {'x':np.array(h5f['dz']['x']), \
             'y':np.array(h5f['dz']['sigma']['y']), \
             'time':np.array(h5f['dz']['sigma']['t']), \
             'z':np.array(h5f['dz']['sigma'][field])})
    D['data']=pc.data().from_h5('/home/ben/temp/'+file, group='data')
    D['data']=D['data'][D['data'].three_sigma_edit==1]
    DD.append(D)
    
fig=plt.figure(); 
plt.clf(); ax=fig.subplots(2,2,gridspec_kw=\
       {'left':0.1, 'right':0.6, 'bottom':0.2, 'top':0.9, 'hspace':0.2, 'wspace':0.2})
ax2=plt.axes(position=[0.65, 0.2, 0.25, 0.7])

colors=['b','g']
for ii in [0, 1]:
    xx=DD[ii]['dz'].x
    yy=np.mean(DD[ii]['dz'].y)+np.zeros_like(xx)
    dzdt=DD[ii]['dzdt_lag2'].interp(xx, yy, band=3)
    sigma_dzdt=DD[ii]['sigma_dzdt_lag2'].interp(xx, yy, band=3)
    plt.axes(ax[ii, 0])
    hi1=DD[ii]['dzdt_lag2'].show(band=3, xy_scale=0.001, cmap='Spectral', vmin=-4, vmax=4)
    plt.plot(D['data'].x/1000, D['data'].y/1000,'.', color='k', markersize=0.25)
    plt.plot(xx/1000, yy/1000,'w')
    plt.axes(ax[ii, 1])
    hi2=DD[ii]['sigma_dzdt_lag2'].show(band=3, xy_scale=0.001, vmin=0.1, vmax=0.25)
    plt.plot(xx/1000, yy/1000,'w')

    plt.axes(ax2)
    plt.errorbar(xx/1000, dzdt, yerr=sigma_dzdt, color=colors[ii])

ax2.yaxis.set_ticks_position('right')
ax2.yaxis.set_label_position("right")
ax2.set_ylabel('elevation rate, m/yr')
for ax_i in ax.flat:
    ax_i.label_outer()


hb1=fig.colorbar(hi1, ax=[ax[0,0], ax[1,0]], orientation='horizontal', fraction=0.1)
hb1.set_label('dz/dt, m/yr')
hb2=fig.colorbar(hi2, ax=[ax[0,1], ax[1,1]], orientation='horizontal', fraction=0.1)
hb2.set_label('dz/dt error, m/yr')
