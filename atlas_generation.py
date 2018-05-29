#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 07:11:16 2018

@author: therry

Module contenant la série d'outils pour les tâches maitre/esclaves telles que :
    - Génération d'atlas
        - Interpolation à l'aide des outils d'interpolation (interpolation_tools.py)
        - Génération des statistiques (zmean, zstd, zdz, ...)
        - Propagation inverse des flag après traitement des données Argo
        - Collage des "dalles" de l'atlas

Ce module contient également les outils permettant la génération d'atlas sur
un domaine restreint :
    - Génération sous-domaine d'atlas avec 4 points du globe (longitude, latitude)
    - Génération sous-domaine dans une zone géographique définie ( Mer Méditerranée,
      Pacifique-Sud, ...)



"""

import param as param
import tile as tile

# Interpolation

def interpolation_on_tiles(itile):
    """
    :param itile: Numéro de la dalle sur laquelle on réalise l'interpolation
    
    Interpolation des profiles ARGO appartenant à la dalle (itile) sur les niveaux
    de référence zref.
    
    Les résultats de l'interpolation sont sauvés dans un fichier pickle.    
    
    :rtype: DataFrame
    """
    return interp_profiles

def compute_stats_at_zref(atlas_infos, itile=None, coord=None):
    """
    :param atlas_infos: Dictionnaire contenant les champs suivants : (mode, date, reso, timeflag, typestat)
                      permettant de connaitre les informations sur les statistiques à calculer
    :param itile: Numéro de la dalle sur laquelle on calcule les statistiques
    :param coord: Dictionnaire contenant les champs suivants : ((lon, lat), deltalon, deltalat)
                  définissant la zone géographique sur laquelle on souhaite calculer les statistiques
    
    Calcul des statistiques sur un domaine défini. Le type de statistique calculé
    dépend du dictionnaire atlas_info passé en argument qui délivre toute l'information
    à la génération de celles-ci.
    Si itile et coord ne sont pas renseignés, la fonction retournera un message
    d'avertissement et produira un résultat sur une dalle géographique prédéfinie.
    La valeur de retour de cette fonction est un dictionnaire contenant les 
    variables statistiques calculées.
    
    :rtype: dict
    """
    stats = {}
    if type(itile) == int and 0 <= itile < 300:
        grid_lat, grid_lon = tile.grid_coordinate(atlas_infos['RESO'], itile = itile)

    elif len(coord) == 3 and len(coord[0]) == 2:
        # on définit les bornes de notre domaine d'étude
        # avec deux latitudes et longitudes (min et max)
        lat = [coord[0][1], (coord[0][1] + coord[2])]
        lon = [coord[0][0], (coord[0][0] + coord[1])]
        points = [lat, lon]
        grid_lat, grid_lon = tile.grid_coordinate(atlas_infos['RESO'], points = points)
    else:
        if itile == None and coord == None:
            itile = 269
            print('No param given (itile or coord)')
            print('Default choice for the evaluation : itile = %i' % itile)
            grid_lat, grid_lon = tile.grid_coordinate(atlas_infos['RESO'], itile)
        else:
            raise TypeError('Bad definition of itile or coord, please refer to the docstring to know how to define them')

    return stats