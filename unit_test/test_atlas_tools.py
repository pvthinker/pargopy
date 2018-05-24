#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:53:00 2018

@author: therry

Module contenant la batterie de tests unitaires pour vérifier le fonctionnement
des fonctions d'atlas_tools.py

"""

import unittest
 
# Le code à tester doit être importable. On
# verra dans une autre partie comment organiser
# son projet pour cela.
from mon_module import get
 
# Cette classe est un groupe de tests. Son nom DOIT commencer
# par 'Test' et la classe DOIT hériter de unittest.TestCase.
class TestFonctionGet(unittest.TestCase):
 
    # Cette méthode sera appelée avant chaque test.
    def setUp(self):
        self.simple_comme_bonjour = ('pomme', 'banane')
 
    # Cette méthode sera appelée après chaque test.
    def tearDown(self):
        print('Nettoyage !')
 
    def test_get_element(self):
        # plus besoin de créer le tuple ici
        element = get(self.simple_comme_bonjour, 0)
        self.assertEqual(element, 'pomme')
 
    def test_element_manquant(self):
        element = get(self.simple_comme_bonjour, 1000, 'Je laisse la main')
        self.assertEqual(element, 'Je laisse la main')
 
    def test_avec_echec(self):
        element = get(self.simple_comme_bonjour, 1000, 'Je laisse la main')
        self.assertEqual(element, 'Je tres clair, Luc')
