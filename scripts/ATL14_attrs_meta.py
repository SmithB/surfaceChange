#!/usr/bin/env python3

import numpy as np
import sys, os
from netCDF4 import Dataset
import h5py
from osgeo import osr, ogr

def write_atl14meta(dst,ncTemplate):

    with Dataset(ncTemplate,'r') as src:
        # copy attributes
        for name in src.ncattrs():
            dst.setncattr(name, src.getncattr(name))
    # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited else None))
    # copy all file data except for the excluded
        for name, variable in src.variables.items():
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst.variables[name][:] = src.variables[name][:]
        for grp in walktree(src):
            for child in grp:
                dg = dst.createGroup(child.path)
                for name in child.ncattrs():
                    dg.setncattr(name,child.getncattr(name))
                for name, dimension in child.dimensions.items():
                    dg.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
                for name, variable in child.variables.items():
                    x = dg.createVariable(name, variable.datatype, variable.dimensions)
                    dg.variables[name][:] = child.variables[name][:]
def walktree(top):
    yield top.groups.values()
    for value in top.groups.values():
        yield from walktree(value)

def set_lineage(dst):
    tilepath = '/att/nobackup/project/icesat-2/ATL14_processing/rel001/north/CN/centers'
# list of lineage attributes
    lineage = []
# regular expression for extracting ATL11 parameters
    rx = re.compile(r'(ATL\d{2})_(\d{4})(\d{2})_(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
# For each tile:
#    for tile in glob.iglob(os.path.join(tilepath,'*.h5')):
# for each file (granule)
#    for FILE in :
    FILE = 'ATL11_000110_0310_004_01.h5'
    if 1 == 1:
# extract parameters from filename
        PRD,TRK,GRAN,SCYC,ECYC,RL,VERS,AUX = rx.findall(FILE).pop()
# extract universally unique identifier from file
        with h5py.File(FILE,'r') as fileID:
            UUID = fileID['METADATA']['DatasetIdentification'].attrs['uuid'].decode('utf-8')
# merge attributes as a tuple
        attrs = (FILE,PRD,int(TRK),int(GRAN),int(SCYC),int(ECYC),int(RL),int(VERS),UUID)
# add attributes to list
        if attrs not in lineage:
            lineage.append(attrs)
# reduce to unique lineage attributes (no repeat files)
    sorted(set(lineage))

def set_geobounds(dst,root_info):
    polar_srs=osr.SpatialReference()
    polar_srs.ImportFromEPSG(int(dst['/Polar_Stereographic'].getncattr('spatial_epsg')))
    ll_srs=osr.SpatialReference()
    ll_srs.ImportFromEPSG(4326)
    if hasattr(osr,'OAMS_TRADITIONAL_GIS_ORDER'):
        ll_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        polar_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ct=osr.CoordinateTransformation(polar_srs, ll_srs)

    xmin,xmax = (np.min(dst['/x']),np.max(dst['/x']))
    ymin,ymax = (np.min(dst['/y']),np.max(dst['/y']))
    N = 2

    multipoint = ogr.Geometryy(ogr.wkbMultiPoint)
    for x in range(N+1):
        for y in range(N+1):
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(ymax - y*dy,xmin + x*dx)
            multipoint.AddGeometry(point)

    multipoint.Transform(ct)
    lonmin,lonmax,latmin,latmax = multipoint.GetEnvelope()
    if (lonmin == -180.0) | (lonmax == 180.0):
        lonmin,lonmax = (-180.0,180.0)
    root_info.update({'geospatial_lon_min': lonmin})
    root_info.update({'geospatial_lon_max': lonmax})
    root_info.update({'geospatial_lat_min': latmin})
    root_info.update({'geospatial_lat_max': latmax})

