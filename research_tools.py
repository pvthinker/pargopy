# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 15:24:20 2018

@author: herry
"""
import os
import numpy as np
import time
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import param as param
import argotools as argotools
import argodb as argodatabase

path_localdata = './data/filter'
argodb = argodatabase.read_argodb()
#  ----------------------------------------------------------------------------
def tile_definition():
    """Creating the variables"""

    minlon = -180.
    maxlon = 180.
    latmin = -80.
    maxlat = 80.

    nlon = 20
    nlat = 15

    minmargin = 1.

    deltalon = maxlon-minlon
    deltalat = maxlat-latmin

    lon = minlon + np.arange(nlon+1) * deltalon/nlon
    lat = latmin + np.arange(nlat+1) * deltalat/nlat

    # distribute the latitudes so that their differences
    # vary in cos(lat)
    # do it via an iterative method
    for k in range(5):
        latm = 0.5*(lat[1:]+lat[:-1])
        dlat = np.diff(lat) * np.cos(latm*np.pi/180)
        dlat = dlat*0 + np.mean(dlat)
        dlat = dlat / np.cos(latm*np.pi/180)
        dlat = dlat/sum(dlat)*deltalat
        lat[1:] = latmin + np.cumsum(dlat)

    margin = minmargin / np.cos(latm*np.pi/180)

    return lat, lon, nlat, nlon, margin


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
                print(k)
                lonmin += 360
            elif lonmax > 180:
                print(k)
                lonmax -= 360
            else:
                pass
            res = {'latmin': latmin,
                   'latmax': latmax,
                   'lonmin': lonmin,
                   'lonmax': lonmax,
                   'lat[i]': lat[i],
                   'lat[i+1]': lat[i+1],
                   'lon[j]': lon[j],
                   'lon[j+1]': lon[j+1],
                   'margin[i]': margin[i]}

            argo_extract = get_idx_from_tiles_lim(res)
            test_tiles(argo_extract, k)
            write_argo_filter(argo_extract, k)
            k += 1


#  ----------------------------------------------------------------------------
def test_tiles(argo_extract, i):
    """ Test to know if the tiles are correctly done with the lat and lon 
    limits"""
    idx1 = np.where(argo_extract['LATITUDE'] > argo_extract['MAXLAT'] + argo_extract['LAT_MARGIN'])
    idx2 = np.where(argo_extract['LATITUDE'] < argo_extract['MINLAT'] - argo_extract['LAT_MARGIN'])
    idx3 = np.where(argo_extract['LONGITUDE'] > argo_extract['MAXLON'] + argo_extract['LON_MARGIN'])
    idx4 = np.where(argo_extract['LONGITUDE'] < argo_extract['MINLON'] - argo_extract['LON_MARGIN'])
    if (idx1[0] != []) | (idx2[0] != []) | (idx3[0] != []) | (idx4[0] != []):
        print('There is an error with the dimensions of the tile number %i' % i)
        exit(0)


#  ----------------------------------------------------------------------------
def write_argo_filter(argo_extract, i):
    with open('%s/argodic%003i.pkl' % (path_localdata, i), 'w') as f:
        pickle.dump(argo_extract, f)


#  ----------------------------------------------------------------------------
def get_idx_from_tiles_lim(res):
    """Get the list of profile indices present in argodb that correspond
       to the list of wmos"""
    #  max and min are the limits with the margins
    if res['lonmin'] > res['lonmax']:
        idx = np.where((argodb['LATITUDE'] > res['latmin']) & (argodb['LATITUDE'] < res['latmax']) & ((argodb['LONGITUDE'] > res['lonmin']) | (argodb['LONGITUDE'] < res['lonmax'])))
    else:
        idx = np.where((argodb['LATITUDE'] > res['latmin']) & (argodb['LATITUDE'] < res['latmax']) & (argodb['LONGITUDE'] > res['lonmin']) & (argodb['LONGITUDE'] < res['lonmax']))
    argo_extract = extract_idx_from_argodb(argodb, idx)
    #  low and high are the limits without margins
    argo_extract['MINLAT'] = res['lat[i]']
    argo_extract['MAXLAT'] = res['lat[i+1]']
    argo_extract['MINLON'] = res['lon[j]']
    argo_extract['MAXLON'] = res['lon[j+1]']
    argo_extract['LAT_MARGIN'] = res['margin[i]']
    argo_extract['LON_MARGIN'] = 2

    return argo_extract


#  ----------------------------------------------------------------------------
def get_idx_from_list_wmo(argodb, wmos):
    """Get the list of profile indices present in argodb that correspond
       to the list of wmos

    """
    infos = argotools.retrieve_infos_from_tag(argodb, argodb['TAG'])
    idx = []
    for w in wmos:
        idx += list(np.where(infos['WMO'] == w)[0])
    return idx


#  ----------------------------------------------------------------------------
def extract_idx_from_argodb(argodb, idx):
    """Return a argodb type dictionnary that is a subset of argodb and
       containing only entries given in idx (list)

    """
    argodb_extract = {}
    for k in argodb.keys():
        argodb_extract[k] = argodb[k][idx]
    return argodb_extract


#  ----------------------------------------------------------------------------
def extract_idx_from_wmostats(wmostats, idx):
    wmostats_extract = {}
    keys = wmostats.keys()
    keys.remove('N_WMO')
    for k in keys:
        wmostats_extract[k] = wmostats[k][idx]
    if type(idx) in [int, np.int64]:
        n_wmo = 1
    else:
        n_wmo = len(idx)
    wmostats_extract['N_WMO'] = n_wmo
    return wmostats_extract


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    creating_tiles()
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
