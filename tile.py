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
import pandas as pd
import pickle as pickle

import param as param
import database as db
import interpolation_tools as interp

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


def retrieve_tile_limits(tiles, itile):
    """
    :param tiles: Dictionnaire contenant les informations relatives à toutes les tiles
    :param itile: Numéro de la tile dont on souhaite connaitre les limites
    
    Fonction retournant les informations pour une tile donnée :
        - lon (min et max)
        - lat (min et max)
        - marginlon
        - marginlat

    :rtype: dict
    """

    j = itile // tiles['NLON']
    i = itile % tiles['NLON']

    lon = [tiles['LONGITUDE'][i], tiles['LONGITUDE'][i+1]]
    lat = [tiles['LATITUDE'][j], tiles['LATITUDE'][j+1]]
    marginlat = tiles['MARGINLAT'][j]
    marginlon = tiles['MARGINLON']

    limits = {'LON_LIM' : lon,
              'LAT_LIM' : lat,
              'MARGINLAT' : marginlat,
              'MARGINLON' : marginlon}

    return limits


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


def read_argo_tile(itile):
    """
    Fonction utilisée pour récupérer la DataFrame 'argo_global' du fichier pickle
    dans lequel elle a été préalablement sauvegardé.

    :rtype: DataFrame
    """
    argo_tile = pd.read_pickle('%s/argo_%003i.pkl' % (param.get_path('parallel'), itile))

    return argo_tile


def write_argo_tile(itile, argo_tile):
    """
    :param itile: Numéro de la tile que l'on veut sauvegarder
    :param argo_tile: Tile à sauvegarder

    Fonction utilisée pour écrire notre DataFrame argo_tile au sein d'un fichier
    pickle portant le nom de : argo_%003i.pkl.

    :rtype: None
    """
    argo_tile.to_pickle('%s/argo_%003i.pkl' % (param.get_path('parallel'), itile))


def write_zref_profiles(itile, zref_profiles):
    """
    :param itile: Numéro de la tile que l'on veut sauvegarder
    :param zref_profiles: Dictionnaire contenant les profiles interpolés

    Fonction utilisée pour écrire nos profiles interpolés dans un fichier pickle

    :rtype: None
    """
    with open('%s/zref_profiles_%003i.pkl' % (param.get_path('zref_profiles'), itile), 'w') as f:
        pickle.dump(zref_profiles, f)


def generate_zref_profiles(itile):
    """Interpolate all Argo profiles in tile 'i' onto 'zref' depths. Save
    the result in the 'tile%003i.pkl' file

    """

    wmodic = db.create_wmodic()
    argo_tile = read_argo_tile(itile)
    zref_profiles = interp.interpolate_profiles(argo_tile, wmodic)
    zref_profiles['ZREF'] = param.zref

    write_zref_profiles(itile, zref_profiles)
    # Mise à jour et sauvegarde d'interp_result dans argo_tile
    interp_result = argo_tile.loc['STATUS'] = True
    synchronize_argo_tile_from_interpolation(itile, argo_tile, interp_result)
    # try to reduce memory leakage when processessing all the tiles
    del(zref_profiles, wmodic)


def synchronize_argo_tile_from_global(itile, argo_tile):
    """
    :param itile: Numéro de la tile que l'on souhaite générer
    :param argo_tile: Version antérieure de l'argo_%3i à synchroniser avec la nouvelle
    Fonction utilisée pour mettre à jour les argo_%3i avec les nouveaux profiles de la base Argo,
    mais également pour le retour d'information après l'interpolation des profiles

    :rtype: DataFrame
    """

    idx_modified_status = []

    tiles = tile_definition()
    limits = retrieve_tile_limits(tiles, itile)
    lat_min = limits['LAT_LIM'][0] - limits['MARGINLAT']
    lat_max = limits['LAT_LIM'][1] + limits['MARGINLAT']
    lon_min = limits['LON_LIM'][0] - limits['MARGINLON']
    lon_max = limits['LON_LIM'][1] + limits['MARGINLON']

    argo_global = db.read_argo_global()
    # Mettre à jour les argo_%3i avec argo_global
    part_of_argo = argo_global[(lat_min <= argo_global['LATITUDE']) &
                               (argo_global['LATITUDE'] <= lat_max) & 
                               (lon_min <= argo_global['LONGITUDE']) &
                               (argo_global['LONGITUDE'] <= lon_max)]
    idx_new_prof = part_of_argo.index.difference(argo_tile.index)
    for tag in argo_tile.index:
        if argo_tile.loc[tag, 'STATUS'] != part_of_argo.loc[tag, 'STATUS']:
            idx_modified_status.append(tag)
    argo_tile.loc[idx_modified_status] = part_of_argo.loc[idx_modified_status]

    if len(idx_new_prof) != 0:
        argo_tile = pd.concat([argo_tile, part_of_argo.loc[idx_new_prof]])

    write_argo_tile(itile, argo_tile)
    print('Tile %003i synchronized from argo_global' % itile)


def synchronize_argo_tile_from_interpolation(itile, argo_tile, interp_result):
    """
    :param itile: Numéro de la tile que l'on souhaite générer
    :param argo_tile: Version antérieure de l'argo_%3i à synchroniser avec la nouvelle
    :param interp_result: Version mise à jour d'argo_tile avec compte rendu interpolation
    Fonction utilisée pour mettre à jour les argo_%3i avec les nouveaux profiles de la base Argo,
    mais également pour le retour d'information après l'interpolation des profiles

    :rtype: DataFrame
    """

    for tag in argo_tile.index:
        if argo_tile.loc[tag, 'FLAG'] != interp_result.loc[tag, 'FLAG']:
            print('Flag of tag %i changed from %i to %i' % (tag, argo_tile.loc[tag, 'FLAG'], interp_result.loc[tag, 'FLAG']))
            argo_tile.loc[tag] = interp_result.loc[tag]
    
    write_argo_tile(itile, argo_tile)
    print('Tile %003i synchronized after interpolation' % itile)


#  ----------------------------------------------------------------------------
if __name__ == '__main__':
    
    for i in range(300):
        argo_tile = read_argo_tile(i)
        synchronize_argo_tile(i, argo_tile, argo_global=True)