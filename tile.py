#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:37:57 2018

@author: therry

Module contenant la série d'outils utilisé pour le découpage de l'atlas en dalles :
    - Découpage en 300 dalles en fonction des longitudes/latitudes
    - 

"""

import numpy as np

def tile_definition():
    """
    Define the tiles coordinates, in the form of a vector of lon and
    lat + their margins

    The tile indexing is

    \|-----+-----+-----+-----+-----|
    \| 280 | 281 | 282 | ... | 299 |
    \|-----+-----+-----+-----+-----|
    \| ... | ... | ... | ... | ... |
    \|-----+-----+-----+-----+-----|
    \|  20 |  21 |  22 | ... |  39 |
    \|-----+-----+-----+-----+-----|
    \|   0 |   1 |   2 | ... |  19 |
    \|-----+-----+-----+-----+-----|

    The returned dictionnary tiles is composed by:
        - lat (vector 1D)
        - lon (vector 1D)
        - nlat ( Quantity of latitude in the tile)
        - nlon ( Quantity of longitude in the tile)
        - marginlat (vector 1D)
        - marginlon (int)

    :rtype: dict

    """

    tiles = {}

    minlon = -180
    maxlon = 180
    minlat = -80
    maxlat = 80

    tiles['NLON'] = 20
    tiles['NLAT'] = 15

    minmargin = 1.

    deltalon = maxlon-minlon
    deltalat = maxlat-minlat

    tiles['LONGITUDE'] = minlon + np.arange(tiles['NLON']+1) * deltalon/tiles['NLON']
    tiles['LATITUDE'] = minlat + np.arange(tiles['NLAT']+1) * deltalat/tiles['NLAT']

    for k in range(5):
        latm = 0.5*(tiles['LATITUDE'][1:]+tiles['LATITUDE'][:-1])
        dlat = np.diff(tiles['LATITUDE']) * np.cos(latm*np.pi/180)
        dlat = dlat*0 + np.mean(dlat)
        dlat = dlat / np.cos(latm*np.pi/180)
        dlat = dlat/sum(dlat)*deltalat
        tiles['LATITUDE'][1:] = minlat + np.cumsum(dlat)

    tiles['MARGINLAT'] = minmargin / np.cos(latm*np.pi/180)
    tiles['MARGINLON'] = 2

    return tiles


def grid_coordinate(reso, itile=None, points=None):
    """ 
    Returns the coordinates of each point of the grid for a given tile or points

    coordinates are round multiples of reso_deg
    reso sets the grid resolution, typically 0.5deg
    
    :rtype: numpy.ndarray, numpy.ndarray""" 

    if itile != None:
        tiles = tile_definition()
        j = itile // tiles['NLON'] # latitude
        i = itile % tiles['NLON'] # longitude
        
        lat = [tiles['LATITUDE'][j], tiles['LATITUDE'][j+1]]
        lon = [tiles['LONGITUDE'][i], tiles['LONGITUDE'][i+1]]
    elif points != None:
        lat = points[0]
        lon = points[1]

    latmin = np.ceil(lat[0]/reso)*reso
    latmax = np.floor(lat[1]/reso)*reso
    lonmin = np.ceil(lon[0]/reso)*reso
    lonmax = np.floor(lon[1]/reso)*reso
    if (lonmax % 1. == 0.):
        # remove rightmost point (it's the leftmost point of the next
        # tile). That's also true for +180 (that is == -180).
        lonmax -= reso

    grid_lat = np.arange(latmin, (latmax+reso), reso)
    grid_lon = np.arange(lonmin, (lonmax+reso), reso)

    return(grid_lat, grid_lon)