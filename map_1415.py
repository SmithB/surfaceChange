#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 20:59:05 2019

@author: ben
"""

import glob
import re
import pointCollection as pc
import matplotlib.pyplot as plt
import numpy as np

thedir='/Volumes/ice2/ben/ATL14_test/IS2/z03xlooser_dt10xlooser_80km'
W_ctr=4.e4


MOG=pc.grid.data().from_geotif('/Volumes/ice1/ben/MOG/2005/mog_2005_1km.tif')

plt.figure()
#MOG.show(vmin=14000, vmax=17000, cmap='gray')
ctr_dict={}

XR=[-440000., -120000.]; YR= [-1840000., -1560000.]

hax=[];
hax.append(plt.gcf().add_subplot(3,2,1))

for count in range(1,7):
    if count>1:
        hax.append(plt.gcf().add_subplot(3,2,count, sharex=hax[0], sharey=hax[0]))
    MOG.show(vmin=14000, vmax=17000, cmap='gray')

for file in glob.glob(thedir+'/*/*.h5'):
    xyc=[int(item)*1.e3 for item in re.compile('E(.*)_N(.*).h5').search(file).groups()]
    if not ((xyc[0]>= XR[0]) and (xyc[0] <= XR[1]) & (xyc[1]>= YR[0]) and (xyc[1] <= YR[1])): 
        continue
    
    temp=pc.grid.data().from_h5( file, group='dz/', field_mapping={'z':'dz'})
    
    xyc=[int(item)*1.e3 for item in re.compile('E(.*)_N(.*).h5').search(file).groups()]
    temp=temp.subset(np.mean(temp.x)+np.array([-W_ctr, W_ctr])/2, np.mean(temp.y)+np.array([-W_ctr, W_ctr])/2)
    for count in range(0, 3):
        hax[count*2].imshow(temp.z[:,:,count+1]-temp.z[:,:,count], extent=temp.extent, cmap='Spectral', vmin=-3, vmax=3, origin='lower')

    temp=pc.grid.data().from_h5( file, group='dz/', field_mapping={'z':'count'})
    temp=temp.subset(np.mean(temp.x)+np.array([-W_ctr, W_ctr])/2, np.mean(temp.y)+np.array([-W_ctr, W_ctr])/2)

    for count in range(0, 3):
        hax[count*2+1].imshow((temp.z[:,:,count+1]>0) & (temp.z[:,:,count]>0), extent=temp.extent, vmin=0, vmax=1, origin='lower')

    ctr_dict[(np.mean(temp.x), np.mean(temp.y))]=file

for count in range(3):
    hax[count*2].set_ylabel(f'$\delta$ z\n {temp.t[count]} to {temp.t[count+1]}')
    hax[count*2+1].set_ylabel(f'coverage\n {temp.t[count]} to {temp.t[count+1]}')
    hax[count*2+1].set_yticks([])

for count in range(2):
    hax[count*2].set_xticks([])


plt.gca().set_xlim(MOG.extent[0:2])
plt.gca().set_ylim(MOG.extent[2:4])




