#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 08:04:37 2018

@author: therry

Module contenant tout les paramètres utiles au fonctionnement du programme :
    - Chemins vers les différents répertoires utilisés :

        - Vers la base de données ARGO
        - Vers le répertoire du projet
        - Vers les fichiers de statistique
        - Vers les fichiers d'atlas
        - ...
    - Variables définissant la/les tâche(s) à réaliser
        - Type de statistique (zmean, zstd, zdz, ...)
        - Domaine d'étude (atlas global, dalles, coordonnées, régions géographiques)
        - Précision de la grille (reso_deg)
        - Date du snapchot choisi pour étudier ARGO
        - Statut des données que l'on souhaite traiter (R, A, D)
        - Filtre (saisonnier, mensuel, ...) sur les données traitées
    - Variables utilisées par tout les modules
        - zref
        - filename (atlas, stats)
        - ...
"""

import numpy as np

#  location is used to know where you are working
#  if you're working on datarmor, location takes the value 'DATARMOR'
#  if you're working on herry-s local, location takes the value 'HERRY'

#  location = 'roullet_lops'
location = 'DATARMOR_TMP'
#  location = 'DATARMOR_FINAL'

def get_atlas_infos():
    """
    Fonction définissant les paramètres de l'atlas que l'on souhaite générer.
    Sa valeur de retour est un dictionnaire contenant les informations suivantes :

        - TYPESTAT
        - RESO
        - DATE
        - MODE
        - TIMEFLAG
    
    :rtype: dict
    """
    atlas_infos = {}
    atlas_infos['TYPESTAT'] = 'zmean'
    atlas_infos['RESO'] = 0.5
    atlas_infos['DATE'] = [2017, 12, 31]
    atlas_infos['MODE'] = 'D'
    atlas_infos['TIMEFLAG'] = 'annual'

    return atlas_infos


zref = np.array([0., 10., 20., 30., 40., 50., 60., 70., 80., 90.,
                 100., 110., 120., 130., 140., 150., 160., 170.,
                 180., 190., 200., 220., 240., 260., 280, 300.,
                 320., 340., 360., 380., 400., 450., 500., 550.,
                 600., 650., 700., 750., 800., 850., 900., 950.,
                 1000., 1050., 1100., 1150., 1200., 1250., 1300.,
                 1350., 1400., 1450., 1500., 1550., 1600., 1650.,
                 1700., 1750., 1800., 1850., 1900., 1950.,
                 2000.])

daclist = ['aoml', 'bodc', 'coriolis', 'csio',
           'csiro', 'incois', 'jma', 'kma',
           'kordi', 'meds', 'nmdis']

def get_path(wanted_path, reso = 0):
    """
    :param wanted_path: Chemin auquel on souhaite accéder
    
    Fonction retournant le chemin auquel on souhaite accéder en fonction du 
    choix défini par 'location' et du 'wanted_path'
    Si on souhaite accéder au répertoire contenant les atlas pour un reso donné,
    on passe la valeur du reso en argument de la fonction en plus du nom 'reso'
    
    :rtype: String
    """
    
    paths = {}

    if location == 'DATARMOR_TMP':
        workdir = '/home2/datawork/therry'
    elif location == 'DATARMOR_FINAL':
        workdir = '/home2/datawork/therry/final/'
    else:
        raise ValueError('This location is not referenced : %s' % location)

    paths['pargopy_output'] = '%s/pargopy_output' % workdir

    #  paths['argo'] = '/datawork/fsi2/coriolis-s/public/co05/co0508/gdac/dac'

    paths['argo'] = '/home/ref-argo/gdac/dac'

    paths['pargopy'] = '/home2/datahome/therry/pargopy/'

    paths['database'] = '%s/database' % paths['pargopy_output']

    paths['parallel'] = '%s/parallel' % paths['database']

    paths['stats'] = '%s/stats_work' % paths['pargopy_output']

    paths['zref_profiles'] = '%s/zref_profiles' % paths['pargopy_output']

    paths['atlas'] = '%s/atlas' % paths['pargopy_output']

    paths['reso'] = '%s/reso_%.2f' % (paths['atlas'], reso)

    if wanted_path in paths:
        return paths[wanted_path]
    else:
        raise ValueError('This path is not referenced : %s' % wanted_path)


def get_atlas_filename():
    """
    Fonction permettant de générer le nom d'un fichier atlas en fonction de ses
    paramètres contenus dans atlas_infos.
    
    :rtype: String
    """

    atlas_infos = get_atlas_infos()

    if False: # many subfolders - short name for the atlas
        atlas_name = '%s_%g_annual' % (atlas_infos['TYPESTAT'], atlas_infos['RESO'])
        ncfile = '%s/%s/%s/%s/%s/%s.nc' % (get_path('atlas'), atlas_infos['RESO'],
                                           atlas_infos['YEAR'], atlas_infos['MODE'],
                                           atlas_infos['TYPESTAT'], atlas_name)
        print(ncfile)

    else: # one folder - long name for the atlas
        ncfile ='%s/%s_%g_%s_%s.nc' % (get_path('atlas'), atlas_infos['TYPESTAT'], 
                                       atlas_infos['RESO'], atlas_infos['YEAR'], 
                                       atlas_infos['MODE'])
        print(ncfile)

    return ncfile


def get_argo_filename(dac, wmo):
    """
    :param dac: DAC du profil recherché
    :param wmo: WMO du profil recherché
    
    Fonction permettant de générer le nom d'un fichier de profiles contenu dans
    ARGO.
    
    :rtype: String
    """
    if type(dac) is int:
        dac = daclist[dac]
    profile_filename = '%s/%s/%i/%i_prof.nc' % (get_path('argo'), dac, wmo, wmo)
    
    return profile_filename
