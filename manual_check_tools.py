#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 07:30:13 2018

@author: therry
"""
import time
import numpy as np

import argotools as tools
import param as param
import interpolation_tools as interpolation
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
def update_stats(itile, tile):
    """Returns the stats without the profile checked as wrong

    :rtype: dict

    """
    reso = 0.5
    mode = 'D'
    date = ['2017', '12', '31']
    stats_mode = ['zmean']
    new_stats = vs.compute_at_zref(itile, reso, mode, date, stats_mode, tile_dict=tile)
    print('new stats generated')

    return new_stats


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    tmps1 = time.time()
    tag = 5904180099
    xlat, xlon = retrieve_coords_from_tag(tag)
    itiles = retrieve_itile_from_coords(xlat, xlon)
    for itile in itiles:
        tile = update_tile(itile, tag)
        stats = update_stats(itile, tile)
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
        