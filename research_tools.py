# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 15:24:20 2018

@author: herry
"""
from netCDF4 import Dataset
import os
import numpy as np
import time
import param as param


#  ----------------------------------------------------------------------------
def variables_creation(nb):
    """Creating the variables"""
    #  Generation of the dimension of the table
    rootgrp = Dataset('/local/tmp/1/herry/filter/argo_tile%003i.nc' % nb, "w", format="NETCDF4")
    rootgrp.createDimension('line_number', None)

    #  Generation of the variables of the table

    update = rootgrp.createVariable("update", "d", ('line_number', ))
    update.long_name = 'last update date'
    update.units = 'datetime'

    juliandays = rootgrp.createVariable("juld", "f4", ('line_number', ))
    juliandays.long_name = 'date'
    juliandays.units = 'julian day'

    prof = rootgrp.createVariable('prof', 'i2', ('line_number', ))
    prof.long_name = 'index of the profile in the profiles dir'
    prof.units = 'dimensionless'

    latitudes = rootgrp.createVariable("lat", "f4", ('line_number', ))
    latitudes.long_name = 'latitude'
    latitudes.units = 'degrees'

    longitudes = rootgrp.createVariable("lon", "f4", ('line_number', ))
    longitudes.long_name = 'longitude'
    longitudes.units = 'degrees'

    dac = rootgrp.createVariable('dac', 'b', ('line_number', ))
    dac.long_name = 'index of the dac'
    dac.units = 'dimensionless'

    wmo = rootgrp.createVariable('wmo', 'i', ('line_number', ))
    wmo.long_name = 'wmoid of the wmo'
    wmo.units = 'dimensionless'

    flag = rootgrp.createVariable('flag', 'i', ('line_number', ))
    flag.long_name = 'flag to know if the profile is OK'
    flag.units = 'dimensionless'

    return (latitudes, longitudes, juliandays, dac, wmo, update, prof, flag, rootgrp)
#  ----------------------------------------------------------------------------
def tile_definition():
    """Creating the variables"""

    minlon = -180.
    maxlon = 180.
    minlat = -80.
    maxlat = 80.

    nlon = 20
    nlat = 15

    minmargin = 1.

    deltalon = maxlon-minlon
    deltalat = maxlat-minlat

    lon = minlon + np.arange(nlon+1) * deltalon/nlon
    lat = minlat + np.arange(nlat+1) * deltalat/nlat

    # distribute the latitudes so that their differences
    # vary in cos(lat)
    # do it via an iterative method
    for k in range(5):
        latm = 0.5*(lat[1:]+lat[:-1])
        dlat = np.diff(lat) * np.cos(latm*np.pi/180)
        dlat = dlat*0 + np.mean(dlat)
        dlat = dlat / np.cos(latm*np.pi/180)
        dlat = dlat/sum(dlat)*deltalat
        lat[1:] = minlat + np.cumsum(dlat)

    margin = minmargin / np.cos(latm*np.pi/180)

    return lat, lon, nlat, nlon, margin

#  ----------------------------------------------------------------------------
def reading_variables(minlat, maxlat, minlon, maxlon):
    """Salvation of the profiles datas"""
    #  Reading the netCDF file fom ARGO
    #  path_to_argotile, path_to_argolog = param.data_research_tools()
    if (os.path.isfile('/local/tmp/1/herry/data/argo_log.nc')):
        f = Dataset('/local/tmp/1/herry/data/argo_log.nc', "r", format="NETCDF4")
        lats = f.variables['lat'][:]
        lons = f.variables['lon'][:]
        flag = f.variables['flag'][:]
        if minlon > maxlon:
            idx = np.where((lats > minlat) & (lats < maxlat) & ((lons > minlon) | (lons < maxlon)) & (flag == 0))
            #  idx = [l for idx, l in enumerate(lons) if (l > minlat) | (l < maxlat) & (lons > minlon) & (lons < maxlon) & (flag == 0)]
        else:
            idx = np.where((lats > minlat) & (lats < maxlat) & (lons > minlon) & (lons < maxlon) & (flag == 0))
            #  idx = [l for idx, l in enumerate(lons) if (l > minlat) & (l < maxlat) & (lons > minlon) & (lons < maxlon) & (flag == 0)]
        latitudes = f.variables['lat'][idx]
        longitudes = f.variables['lon'][idx]
        prof = f.variables['prof'][idx]
        update = f.variables['update'][idx]
        juld = f.variables['juld'][idx]
        dac = f.variables['dac'][idx]
        wmo = f.variables['wmo'][idx]
        flag = f.variables['flag'][idx]

        f.close()
    return (latitudes, longitudes, prof, update, juld, dac, wmo, flag)

#  ----------------------------------------------------------------------------
def writing_variables(minlat, maxlat, minlon, maxlon, i):
    """Giving values to the variables"""
    #  Generation of the dimension of import matplotlib.pyplot as plt
    latitudes, longitudes, juliandays, dac1, wmo1, update1, prof1, flag1, rootgrp = variables_creation(i)
    rootgrp.setncattr('latmin', minlat)
    rootgrp.setncattr('latmax', maxlat)
    rootgrp.setncattr('lonmin', minlon)
    rootgrp.setncattr('lonmax', maxlon)
    lats, lons, prof, update, juld, dac, wmo, flag = reading_variables(minlat, maxlat, minlon, maxlon)
    profile_per_tile(len(prof), i)
    latitudes[:] = lats
    longitudes[:] = lons
    juliandays[:] = juld
    flag1[:] = flag
    dac1[:] = dac
    wmo1[:] = wmo
    update1[:] = update
    prof1[:] = prof
    rootgrp.close()
    print('argo_tile%003i created' % i)


#  ----------------------------------------------------------------------------
def creating_tiles():
    """Giving values to the variables"""
    #  Generation of the dimension of import matplotlib.pyplot as plt
    lat, lon, nlat, nlon, margin = tile_definition()
    k = 0
    for i in range(nlat):
        for j in range(nlon):
            latmin = lat[i] - margin[i]
            latmax = lat[i + 1] + margin[i]
            lonmin = lon[j] - 2
            lonmax = lon[j + 1] + 2
            if lonmin < -180:
                lonmin += 360
            elif lonmax > 180:
                lonmax -= 360
            else:
                pass
            writing_variables(latmin, latmax, lonmin, lonmax, k)
            k += 1

#  ----------------------------------------------------------------------------
def retrieve_dac_from_wmoid(dac, wmo, wmoid):
    """Retrieving the dac with the wmo"""
    i = 0
    while wmo[i] != wmoid:
        i += 1
    return(dac[i])


#  ----------------------------------------------------------------------------
def retrieve_prof_from_wmoid(prof, wmo, wmoid):
    """Retrieving the profiles index corresponding to the same the wmo"""
    prof1 = []
    for i, profile in enumerate(prof):
        if wmo[i] == wmoid:
            prof1.append(prof[i])
    return(prof1)

#  ----------------------------------------------------------------------------
def profile_per_tile(nb_prof, i):
    """Counting the number of profile per tile"""
    profile = np.zeros((300,))
    profile[i] = nb_prof
    print('There are %i profiles in this tile' % profile[i])
    return profile

#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    creating_tiles()
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
