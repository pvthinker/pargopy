#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 07:11:16 2018

@author: therry

Module contenant la série d'outil pour les tâches maitre/esclaves telles que :
    - Génération d'atlas
        -> Interpolation
        -> Génération des statistiques (zmean, zstd, zdz, ...)
        -> Propagation inverse des flag après traitement des données Argo
        -> Collage des "dalles" de l'atlas

Ce module contient également les outils permettant la génération d'atlas sur
un domaine restreint :
    - Génération sous-domaine d'atlas avec 4 points du globe (longitude, latitude)
    - Génération sous-domaine dans une zone géographique définie ( Mer Méditerranée,
      Pacifique-Sud, ...)



"""

import param as param

