#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:53:00 2018

@author: therry

Module contenant la batterie de tests unitaires pour vérifier le fonctionnement
des fonctions de general_tools.py

"""
import os
import sys
sys.path[:0] = ['../']
import tile as tile
import param as param


def test_tile_definition():
    tiles = tile.tile_definition()
    assert len(tiles['LATITUDE']) == tiles['NLAT'] + 1
    assert len(tiles['LONGITUDE']) == tiles['NLON'] + 1


#==============================================================================
# def test_grid_coordinate_with_itile():
#     grid_lat, grid_lon = tile.grid_coordinate(0.5, itile=269)
#     nb_lat = 0 # a renseigner !!!
#     nb_lon = 0 # a renseigner !!!
#     assert len(grid_lat) == nb_lat
#     assert len(grid_lon) == nb_lon
#==============================================================================


def test_grid_coordinate_with_points():
    lat = [10, 20]
    lon = [10, 20]
    reso = 0.5
    grid_lat, grid_lon = tile.grid_coordinate(reso, points=[lat, lon])
    nb_lat = (lat[1] - lat[0])/reso + 1
    # grid_coordinate() remove the last point of grid_lon (cf tile.py)
    nb_lon = (lon[1] - lon[0])/reso
    assert len(grid_lat) == nb_lat
    assert len(grid_lon) == nb_lon

              
def test_zref_profiles():
    not_file = 0
    not_file_list = []
    for itile in range(300):
        if os.path.isfile('%s/zref_profiles_%003i.pkl' % (param.get_path('zref_profiles'), itile)):
            pass
#==============================================================================
#             zref_profile = tile.read_zref_profiles(itile)
#             if len(zref_profile['CT']) == 0:
#                 print('Tile %i empty' % itile)
#==============================================================================
            
        else:
            not_file += 1
            not_file_list.append(itile)
    print('%i file missing' % not_file)
    print('Id of tiles missing', not_file_list)

def test_argo_tile():
    for itile in range(300):
        argo_tile = tile.read_argo_tile(itile)
        if len(argo_tile.index) == 0:
                print('Tile %i empty' % itile)