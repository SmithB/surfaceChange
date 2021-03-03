#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 21:00:02 2019

@author: ben
"""

#import matplotlib.pyplot as plt
import numpy as np
import pointCollection as pc
import scipy.ndimage as snd
import sys
import os
import re
import argparse

def pad_mask_canvas(D, N):
    dx=np.diff(D.x[0:2])
    left=np.arange(-N*dx,0, dx)
    right=np.arange(0, N*dx, dx)
    x1=np.unique(np.concatenate([left+D.x[0], D.x, D.x[-1]+right]))
    y1=np.unique(np.concatenate([left+D.y[0], D.y, D.y[-1]+right]))
    cols=np.flatnonzero(np.in1d(x1, D.x))
    rows=np.flatnonzero(np.in1d(y1, D.y))
    z1=np.zeros([y1.size, x1.size], dtype='bool')
    z1[rows[0]:rows[-1]+1,cols[0]:cols[-1]+1]=D.z.astype('bool')
    return pc.grid.data().from_dict({'x':x1, 'y':y1,'z':z1})

# define the script
prog = "/home/besmith4/git_repos/surfaceChange/surfaceChange/ATL11_to_ATL15.py"

# account for a bug in argparse that misinterprets negative agruents
argv=sys.argv
for i, arg in enumerate(argv):
    if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg


parser=argparse.ArgumentParser(description="generate a list of commands to run ATL11_to_ATL15")
parser.add_argument('step', type=str)
parser.add_argument('defaults_file', type=str)
parser.add_argument('--region_file', '-R', type=str)
parser.add_argument('--skip_errors','-s', action='store_true')
args=parser.parse_args()

if args.step not in ['centers', 'edges','corners']:
    raise(ValueError('argument not known'))
    sys.exit()

if args.skip_errors:
    calc_errors=False
else:
    calc_errors=True
    
XR=None
YR=None
if args.region_file is not None:
    line_re=re.compile('(..)\s*=\s*\[\s*(\S+),\s*(\S+)\s*]')
    temp={}
    with open(args.region_file,'r') as fh:
        for line in fh:
            m = line_re.search(line)
            temp[m.group(1)]=[np.float(m.group(2)), np.float(m.group(3))]
    XR=temp['XR']
    YR=temp['YR']

defaults_re=re.compile('(.*)\s*=\s*(.*)')
defaults={}
with open(args.defaults_file,'r') as fh:
   for line in fh:
       m=defaults_re.search(line)
       if m is not None:
           defaults[m.group(1)]=m.group(2)

try:
    out_dir=defaults['-b']
except keyError:
    out_dir=defaults['--base_dir']

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
step_dir=out_dir+'/'+args.step
if not os.path.isdir(step_dir):
    os.mkdir(step_dir)


Wxy=float(defaults['-W'])
Hxy=Wxy/2

mask_base, mask_ext = os.path.splitext(defaults['--mask_file'])
if mask_ext=='.tif': 
    temp=pc.grid.data().from_geotif(defaults['--mask_file'].replace('100m','1km'))
    mask_G=pad_mask_canvas(temp, 200)
    mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([1, int(3*Hxy/1000)+1], dtype='bool'))
    mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([int(3*Hxy/1000)+1, 1], dtype='bool'))
    x0=np.unique(np.round(mask_G.x/Hxy)*Hxy)
    y0=np.unique(np.round(mask_G.y/Hxy)*Hxy)
    x0, y0 = np.meshgrid(x0, y0)
    xg=x0.ravel()
    yg=y0.ravel()
    good=(np.abs(mask_G.interp(xg, yg)-1)<0.1) & (np.mod(xg, Wxy)==0) & (np.mod(yg, Wxy)==0)

elif mask_ext in ['.shp','.db']:
    mask_G=pc.grid.data().from_geotif(mask_base+'_80km.tif')
    xg, yg = np.meshgrid(mask_G.x, mask_G.y)
    xg=xg.ravel()[mask_G.z.ravel()==1]
    yg=yg.ravel()[mask_G.z.ravel()==1]
    good=np.ones_like(xg, dtype=bool)

if XR is not None:
    good &= (xg>=XR[0]) & (xg <= XR[1]) & (yg > YR[0]) & (yg < YR[1])

xg=xg[good]
yg=yg[good]

    
if args.step=='centers':
    delta_x=[0]
    delta_y=[0]
    symbol='r*'
elif args.step=='edges':
    delta_x=[-1, 0, 0, 1.]
    delta_y=[0, -1, 1, 0.]
    symbol='bo'
elif args.step=='corners':
    delta_x=[-1, 1, -1, 1.]
    delta_y=[-1, -1, 1, 1.]
    symbol='ms'

queued=[];
for xy0 in zip(xg, yg):
    for dx, dy in zip(delta_x, delta_y):  
        xy1=np.array(xy0)+np.array([dx, dy])*Hxy
        if tuple(xy1) in queued:
            continue
        else:
            queued.append(tuple(xy1))
        out_file=out_dir+'/%s//E%d_N%d.h5' % (args.step, xy1[0]/1000, xy1[1]/1000)  
        if not os.path.isfile(out_file):
            #plt.plot(xy0[0], xy0[1],'r*')
            #if np.sqrt((xy1[0]-320000.)**2 + (xy1[1]- -2520000.)**2) > 2.e5:
            #    continue
            #plt.plot(xy1[0], xy1[1],symbol)
            cmd='python3 %s %d %d --%s @%s ' % (prog, xy1[0], xy1[1], args.step, args.defaults_file)
            if calc_errors:
                cmd += ';'+cmd+' --calc_error_for_xy'
            print('source activate IS2; '+cmd)
            

