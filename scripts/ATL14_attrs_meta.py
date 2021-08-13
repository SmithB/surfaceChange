#!/usr/bin/env python3

import numpy as np
import sys, os
from netCDF4 import Dataset
import netCDF4
import h5py
from osgeo import osr, ogr
import csv
import re
import glob
import uuid
from datetime import datetime

def write_atl14meta(dst,fileout,ncTemplate):

    root_info={'asas_release':'SET_BY_PGE', 'date_created':'', 'fileName':'', 'geospatial_lat_max':0., \
        'geospatial_lat_min':0., 'geospatial_lon_max':0., 'geospatial_lon_min':0., \
        'netcdfversion':'', 'history':'SET_BY_PGE', \
        'identifier_product_format_version':'SET_BY_PGE', 'time_coverage_duration':0., \
        'time_coverage_end':'', 'time_coverage_start':'', 'uuid':''}
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
    set_lineage(dst,root_info)
    set_geobounds(dst,root_info)

# to get netcdf version netCDF4.__netcdf4libversion__
    root_info.update({'netcdfversion': netCDF4.__netcdf4libversion__})
    root_info.update({'uuid': str(uuid.uuid4())})
    dst['METADATA/DatasetIdentification'].setncattr('uuid', str(uuid.uuid4()).encode('ASCII'))
    root_info.update({'date_created': str(datetime.now().date())})
    dst['METADATA/DatasetIdentification'].setncattr('creationDate', str(datetime.now().date()))
    root_info.update({'fileName': os.path.basename(fileout)})
    dst['METADATA/DatasetIdentification'].setncattr('fileName', os.path.basename(fileout))
    print(root_info)
    for key, keyval in root_info.items():
        dst.setncattr(key, keyval)

def copy_group(dst,child):
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

def set_lineage(dst,root_info):
    tilepath = '/att/nobackup/project/icesat-2/ATL14_processing/rel001/north/CN/centers'
    atl11path = '/att/nobackup/project/icesat-2/ATL14_processing/ATL11_rel004/north'
    firstATL11 = True
# list of lineage attributes
    lineage = []
# regular expression for extracting ATL11 parameters
    rx = re.compile(r'(ATL\d{2})_(\d{4})(\d{2})_(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
# For each tile:
    min_start_delta_time = np.finfo(np.float64()).max
    max_end_delta_time = np.finfo(np.float64()).tiny
    for tile in glob.iglob(os.path.join(tilepath,'*.h5')):
        with h5py.File(tile,'r') as cf:
            
# for each file (granule)
          for FILE in cf['/meta'].attrs['input_files'].split(','):
# extract parameters from filename
            PRD,TRK,GRAN,SCYC,ECYC,RL,VERS,AUX = rx.findall(FILE).pop()
# extract universally unique identifier from file
            with h5py.File(os.path.join(atl11path,FILE),'r') as fileID:
                UUID = fileID['METADATA']['DatasetIdentification'].attrs['uuid'].decode('utf-8')
                SGEOSEG = fileID['ancillary_data/start_geoseg'][0]
                EGEOSEG = fileID['ancillary_data/end_geoseg'][0]
                SORBIT = fileID['ancillary_data/start_orbit'][0]
                EORBIT = fileID['ancillary_data/end_orbit'][0]
                sdeltatime = fileID['ancillary_data/start_delta_time'][0]
                edeltatime = fileID['ancillary_data/end_delta_time'][0]
                if sdeltatime < min_start_delta_time:
                    sUTCtime = fileID['ancillary_data/data_start_utc'][0].decode('utf-8')
                    min_start_delta_time = sdeltatime
                if edeltatime > max_end_delta_time:
                    eUTCtime = fileID['ancillary_data/data_end_utc'][0].decode('utf-8')
                    max_end_delta_time = edeltatime

                
# merge attributes as a tuple
            attrs = (FILE,PRD,int(TRK),int(GRAN),int(SCYC),int(ECYC),int(VERS),UUID,int(SGEOSEG),
                    int(EGEOSEG),int(SORBIT),int(EORBIT))
# add attributes to list
            if attrs not in lineage:
                lineage.append(attrs)
# reduce to unique lineage attributes (no repeat files)
#    sorted(set(lineage))
    slineage = sorted(lineage,key=lambda x: (x[0]))
    dst['METADATA/Lineage/ATL11'].setncattr('fileName',list(zip(*slineage))[0])
    dst['METADATA/Lineage/ATL11'].setncattr('shortName',list(zip(*slineage))[1])
    dst['METADATA/Lineage/ATL11'].setncattr('start_rgt',list(zip(*slineage))[2])
    dst['METADATA/Lineage/ATL11'].setncattr('end_rgt',list(zip(*slineage))[2])
    dst['METADATA/Lineage/ATL11'].setncattr('start_region',list(zip(*slineage))[3])
    dst['METADATA/Lineage/ATL11'].setncattr('end_region',list(zip(*slineage))[3])
    dst['METADATA/Lineage/ATL11'].setncattr('start_cycle',list(zip(*slineage))[4])
    dst['METADATA/Lineage/ATL11'].setncattr('end_cycle',list(zip(*slineage))[5])
    dst['METADATA/Lineage/ATL11'].setncattr('version',list(zip(*slineage))[6])
    dst['METADATA/Lineage/ATL11'].setncattr('uuid',list(zip(*slineage))[7])
    dst['METADATA/Lineage/ATL11'].setncattr('start_geoseg',list(zip(*slineage))[8])
    dst['METADATA/Lineage/ATL11'].setncattr('end_geoseg',list(zip(*slineage))[9])
    dst['METADATA/Lineage/ATL11'].setncattr('start_orbit',list(zip(*slineage))[10])
    dst['METADATA/Lineage/ATL11'].setncattr('end_orbit',list(zip(*slineage))[11])

    print('start/end/delta time',sUTCtime,eUTCtime,max_end_delta_time-min_start_delta_time)
    root_info.update({'time_coverage_start': sUTCtime})
    root_info.update({'time_coverage_end': eUTCtime})
    root_info.update({'time_coverage_duration': max_end_delta_time-min_start_delta_time})
    dst['/METADATA/Extent'].setncattr('rangeBeginningDateTime',sUTCtime)
    dst['/METADATA/Extent'].setncattr('rangeEndingDateTime',eUTCtime)

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
    dx = (xmax-xmin)/N
    dy = (ymax-ymin)/N

    multipoint = ogr.Geometry(ogr.wkbMultiPoint)
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

    print('lonmin,lonmax,latmin,latmax:',lonmin,lonmax,latmin,latmax)
#    if '/orbit_info/bounding_polygon_dim1' not in dst:
#        dst.createVariable('/orbit_info/bounding_polygon_dim1',

#    dst.variables['/orbit_info/bounding_polygon_dim1'][:] = np.arange(1,4)
#    lon1 = dst['/orbit_info'].variables['bounding_polygon_lon1'][:]
#    print('lon1:',lon1)

    dst['/orbit_info'].variables['bounding_polygon_dim1'][:] = np.arange(1,4+1)
    dst['/orbit_info'].variables['bounding_polygon_lon1'][:] = np.array([lonmin,lonmax,lonmax,lonmin])[:]
    dst['/orbit_info'].variables['bounding_polygon_lat1'][:] = np.array([latmax,latmax,latmin,latmin])[:]
    dst['/METADATA/Extent'].setncattr('westBoundLongitude',lonmin)
    dst['/METADATA/Extent'].setncattr('eastBoundLongitude',lonmax)
    dst['/METADATA/Extent'].setncattr('northBoundLatitude',latmax)
    dst['/METADATA/Extent'].setncattr('southBoundLatitude',latmin)

