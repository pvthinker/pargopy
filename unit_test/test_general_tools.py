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
import general_tools as tools

def test_conversion_gregd_juld():
    gregday = tools.conversion_juld_gregd(19000.5)
    julday = tools.conversion_gregd_juld(gregday)
    assert julday == 19000.5
 
def test_tile_definition():
    tiles = tools.tile_definition()
    assert len(tiles['LATITUDE']) == tiles['NLAT']
    assert len(tiles['LONGITUDE']) == tiles['NLON']
