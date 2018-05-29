#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:53:00 2018

@author: therry

Module contenant la batterie de tests unitaires pour vérifier le fonctionnement
des fonctions de general_tools.py

"""
import sys
sys.path[:0] = ['../']
import tile as tile


def test_tile_definition():
    tiles = tile.tile_definition()
    assert len(tiles['LATITUDE']) == tiles['NLAT']
    assert len(tiles['LONGITUDE']) == tiles['NLON']


def test_grid_coordinate_with_itile():
    grid_lat, grid_lon = tile.grid_coordinate(0.5, itile=269)
    nb_lat = 0 # a renseigner !!!
    nb_lon = 0 # a renseigner !!!
    assert len(grid_lat) == nb_lat
    assert len(grid_lon) == nb_lon


def test_grid_coordinate_with_points():
    lat = [10, 20]
    lon = [10, 20]
    reso = 0.5
    grid_lat, grid_lon = tile.grid_coordinate(reso, points=[lat, lon])
    nb_lat = (lat[1] - lat[0])/reso + 1
    nb_lon = (lon[1] - lon[0])/reso + 1
    assert len(grid_lat) == nb_lat
    assert len(grid_lon) == nb_lon