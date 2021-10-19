import time
import matplotlib.pyplot as plt
import numpy as np
import glob
import h5py
import os
import re
import pointCollection as pc

class tile_picker(object):
    def __init__(self, thedir, handles=None, W=8.e4, map_data=None, **map_args):
        
        tile_re=re.compile('E(.*)_N(.*).h5')
        
        self.xy_file_dict = {tuple(1000*np.array([*map(int, tile_re.search(ff).groups())])):ff 
                             for ff in glob.glob(thedir+'/*/E*_N*.h5') }
        self.xy_tiles = np.array(list(self.xy_file_dict.keys()))
        print(self.xy_tiles.shape)

        if handles is not None and len(handles):
            self.handles=handles
        else:
            self.handles={}
            self.__init_new_ui__(map_data, map_args)
        self.messages=[[]]
        self.last_pt=[[]]
        self.last_file=''
        self.dz_picker=None
        self.last_click_time=0.0
        self.max_click_time = 0.1

        self.W=W

        self.cid=self.handles['figure'].canvas.mpl_connect('button_press_event', self.buttondown)
        self.cid=self.handles['figure'].canvas.mpl_connect('button_release_event', self.buttonup)
    
    def __init_new_ui__(self, map_data, map_args):
        if 'figure' not in self.handles:
            self.handles['figure']=plt.figure()
        if 'tiles_ax' not in self.handles:
            self.handles['tiles_ax'], self.handles['messages']=self.handles['figure'].subplots(1,2)
        if map_data is not None:
            map_data.show(ax=self.handles['tiles_ax'], **map_args)
        self.handles['tiles_ax'].plot(self.xy_tiles[:,0], self.xy_tiles[:,1],'k.')
            
    def buttondown(self, event):
        if not event.inaxes in [self.handles['tiles_ax']]:
            return
        self.last_click_time=time.time()
    def buttonup(self, event):
        try:
            if not event.inaxes in [self.handles['tiles_ax']]:
                self.messages += ['tile_picker: last point not in tiles axis']
                return
            dt_click = time.time()-self.last_click_time
            if time.time()-self.last_click_time > self.max_click_time:
                self.messages += [f'too much time has elapsed : {dt_click}']
                return
            xy0=(event.xdata, event.ydata)
            xy_tile = tuple((np.round(np.array(xy0)/(self.W/2))*self.W/2).astype(int))
            self.messages = [f'xy0={xy0}, xy_tile={xy_tile}']
            if xy_tile not in self.xy_file_dict:
                self.messages += [f'searching by dist for {xy0}']
                this = np.argmin((self.xy_tiles[:,0]-xy0[0])**2 + (self.xy_tiles[:,1]-xy0[1])**2)
                xy_tile = tuple(self.xy_tiles[this,:]) 
            self.last_file=self.xy_file_dict[xy_tile]
            self.handles['tiles_ax'].plot(xy0[0], xy0[1],'x')
            self.handles['tiles_ax'].plot(xy_tile[0], xy_tile[1],'r.')

        except Exception as e:
            self.messages += [e]
            self.handles['tiles_ax'].set_title('ERROR')