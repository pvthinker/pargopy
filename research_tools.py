# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 15:24:20 2018

@author: herry
"""
import os
import numpy as np
import time
import pickle
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
            idx = np.where((argodb['LATITUDE'] > latmin) & (argodb['LATITUDE'] < latmax) & (argodb['LONGITUDE'] > lonmin) & (argodb['LONGITUDE'] < lonmax))
            argo_extract = extract_idx_from_argodb(argodb, idx)
            argo_extract['MINLAT'] = lat[i]
            argo_extract['MAXLAT'] = lat[i + 1]
            argo_extract['MINLON'] = lon[i]
            argo_extract['MAXLON'] = lon[i + 1]
            argo_extract['MARGIN'] = margin[i]
            write_argo_filter(argo_extract, k)
            k += 1


#  ----------------------------------------------------------------------------
def write_argo_filter(argo_extract, i):
    with open('%s/argodic%003i.pkl' % (path_localdata, i), 'w') as f:
        pickle.dump(argo_extract, f)


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
