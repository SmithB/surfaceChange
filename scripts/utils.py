#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 11:12:33 2020

@author: suzanne
"""

def projection_details(fhandle, hemisphere):
    if hemisphere == 1:
        fhandle.attrs['Geographic Coordinate System'] = 'WGS 84'
        fhandle.attrs['Projected Coordinate System'] = 'WGS 84 / NSIDC Sea Ice Polar Stereographic North'
        fhandle.attrs['Longitude of True Origin'] = '-45'
        fhandle.attrs['Latitude of True Origin'] = '70'
        fhandle.attrs['Scale factor at longitude of true origin'] = '1'
        fhandle.attrs['Datum'] = 'WGS 1984'
        fhandle.attrs['Ellipsoid/spheroid'] = 'WGS 84'
        fhandle.attrs['Units'] = 'meter'
        fhandle.attrs['False Easting'] = '0'
        fhandle.attrs['False Northing'] = '0'
        fhandle.attrs['PROJ4 String'] = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
        fhandle.attrs['EPSG Code'] = 'http://epsg.io/3413'
    else:
        fhandle.attrs['Geographic Coordinate System'] = 'WGS 84'
        fhandle.attrs['Projected Coordinate System'] = 'WGS 84 / Antarctic Polar Stereographic'            
        fhandle.attrs['Longitude of True Origin'] = '0'
        fhandle.attrs['Latitude of True Origin'] = '-71'
        fhandle.attrs['Scale factor at longitude of true origin'] = '1'
        fhandle.attrs['Datum'] = 'WGS 1984'
        fhandle.attrs['Ellipsoid/spheroid'] = 'WGS 84'
        fhandle.attrs['Units'] = 'meter'
        fhandle.attrs['False Easting'] = '0'
        fhandle.attrs['False Northing'] = '0'
        fhandle.attrs['PROJ4 String'] = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
        fhandle.attrs['EPSG Code'] = 'http://epsg.io/3031'
    return fhandle
    