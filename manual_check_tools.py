#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 07:30:13 2018

@author: therry
"""
import time
import numpy as np
from netCDF4  import Dataset

import argotools as tools
import param as param
import interpolation_tools as interpolation
import stats as stats
import variable_selector as vs

path_to_data = param.path_to_data
path_to_filter = param.path_to_filter
path_to_tiles = param.path_to_tiles

zref = tools.zref

#  ----------------------------------------------------------------------------
def retrieve_coords_from_tag(tag):
    """Retrieve the coords (lon, lat) of a given tag of a profile

    :rtype: list[int, int]
    """
    argodb = tools.read_dic('argodb', path_to_data)
    idx = np.where(argodb['TAG'] == tag)
    xlat = argodb['LATITUDE'][idx][0]
    xlon = argodb['LONGITUDE'][idx][0]
    
    return xlat, xlon

#  ----------------------------------------------------------------------------
def retrieve_itile_from_coords(xlat, xlon):
    """Retrieve the tile which contains the profile with its lat and lon

    :rtype: int

    """
    list_of_itiles = []
    lat, lon, nlat, nlon, marginlat, marginlon = tools.tile_definition()
    for i in range(nlat):
        for j in range(nlon):
            itile = 20 * i + j
            if lat[i] < xlat and xlat < lat[i+1]:
                if lon[j] < xlon and xlon < lon[j+1]:
                    print('This profile corresponds to tile%003i' % itile)
                    list_of_itiles.append(itile)
    return list_of_itiles


#  ----------------------------------------------------------------------------
def update_argoextract(itile, tag):
    """Returns an argo_extract dict where the profile flagged as wrong with the
    given tag is updated (flag updated with value 404)

    :rtype: dict

    """
    argoextract = tools.read_dic('argodic%003i' % itile, path_to_filter)
    idx = np.where(argoextract['TAG'] == tag)
    print('old argoextract flag : %i' % argoextract['FLAG'][idx[0]])
    argoextract['FLAG'][idx[0]] = '404'
    print('new argoextract flag : %i' % argoextract['FLAG'][idx[0]])
    return argoextract


#  ----------------------------------------------------------------------------
def update_tile(itile, tag):
    """Returns the tile updated without the tag checked as wrong

    :rtype: dict

    """
    new_tile = {}
    keys = ['JULD', 'LONGITUDE', 'DATA_MODE', 'TAG', 'BVF2', 'RHO', 'LATITUDE', 'SA', 'CT']
    tile = tools.read_dic('tile%003i' % itile, path_to_tiles)
    idx = np.where(tile['TAG'] != tag)
    for k in keys:
        new_tile[k] = tile[k][idx]
    new_tile['ZREF'] = zref
    print('new tile generated')

    return new_tile


#  ----------------------------------------------------------------------------
def calculate_new_stats(itile, tile, xlat, xlon):
    """This function returns a dictionnary of all the stats calculated on only
    9 points around the xlat, xlon given in parameters.
    Its ame is to calculate the bew stats without the given bad profile known
    with his coordinates (xlat, xlon).
    X     X     X
    
        +
    X     X     X
    
    
    X     X     X
    
    + : the given coordinate which we found the bad profile
    X : the points that will be calculated after deleting the bad profile
    :rtype: dict

    """
    reso = 0.5
    mode = 'D'
    date = ['2017', '12', '31']
    #  stats_mode = ['zmean']
    #  new_stats = vs.compute_at_zref(itile, reso, mode, date, stats_mode, tile_dict=tile)
    grid_lat, grid_lon = stats.grid_coordinate(itile, reso)
    ilat = int(xlat/reso)*reso
    ilon = int(xlon/reso)*reso
    idx_lat = np.where(grid_lat == ilat)
    idx_lon = np.where(grid_lon == ilon)
    near_lats = []
    near_lons = []
    nb = [-1, 0, 1]
    for i in nb:
        near_lats.append(grid_lat[idx_lat[0][0]+i])
        near_lons.append(grid_lon[idx_lon[0][0]+i])
    print(near_lats, near_lons)
    new_stats = stats.compute_stats_at_zref(mode, date, near_lons, near_lats, reso, manual_check=True, new_tile=tile)
    print('new stats generated')
    return new_stats


#  ----------------------------------------------------------------------------
def update_stats(var, var_name, new_stats):
    """Returns the given variable updated without the value which has the 
    tag checked as wrong

    :rtype: dict

    """
    reso = 0.5
    itile = 1
    idx_j = []
    idx_i = []
    grid_lat, grid_lon = stats.grid_coordinate(itile, reso)
    lon_deg, lat_deg = np.meshgrid(grid_lon, grid_lat)
    nlat, nlon = np.shape(lon_deg)
    for k in range(len(new_stats['lat'])):
        idx_lat = np.where(lat_deg == new_stats['lat'][k])
        idx_lon = np.where(lon_deg[0] == new_stats['lon'][k])
        idx_j.append(idx_lat[0][0])
        idx_i.append(idx_lon[0][0])
    
    for j, idxj in enumerate(idx_j):
        for i, idxi in enumerate(idx_i):
            print(var[:, idxj, idxi])
            var[:, idxj, idxi] = new_stats[var_name][:, j, i]
            print('New stats :')
            print(var[:, idxj, idxi])

    return var
    
#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    
    ncfile = '/home2/datawork/therry/tmp/atlas/0.5/2017/D/zstd/zstd_0.5_annual.nc'

    with Dataset(ncfile, 'r', format='NETCDF4') as nc:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        # zref = nc.variables['zref'][:]
        var =  nc.variables['SAstd'][:, :, :]
    
    tag = 5904180099
    xlat, xlon = retrieve_coords_from_tag(tag)
    itiles = retrieve_itile_from_coords(xlat, xlon)
    for itile in itiles:
        tile = update_tile(itile, tag)
        new_stats = calculate_new_stats(itile, tile, xlat, xlon)
        new_var = update_stats(var, 'SAstd', new_stats)
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
    

        