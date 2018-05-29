#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:53:00 2018

@author: therry

Module contenant la batterie de tests unitaires pour vérifier le fonctionnement
des fonctions d'atlas_generation.py

"""
import sys
sys.path[:0] = ['../']
import atlas_generation as at_gen
 
def test_get():
    simple_comme_bonjour = ('pomme', 'banane')
    element = get(simple_comme_bonjour, 0)
    assert element == 'pomme'
 
def test_element_manquant():
    simple_comme_bonjour = ('pomme', 'banane')
    element = get(simple_comme_bonjour, 1000, 'Je laisse la main')
    assert element == 'Je laisse la main'