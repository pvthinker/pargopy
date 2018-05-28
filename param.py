#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 08:04:37 2018

@author: therry

Module contenant tout les paramètres utiles au fonctionnement du programme :
    - Chemins vers les différents répertoires utilisés
        -> Vers la base de données ARGO
        -> Vers le répertoire du projet
        -> Vers les fichiers de statistique
        -> Vers les fichiers d'atlas
        -> ...
    - Variables définissant la/les tâche(s) à réaliser
        -> Type de statistique (zmean, zstd, zdz, ...)
        -> Domaine d'étude (atlas global, dalles, coordonnées, régions géographiques)
        -> Précision de la grille (reso_deg)
        -> Date du snapchot choisi pour étudier ARGO
        -> Statut des données que l'on souhaite traiter (R, A, D)
        -> Filtre (saisonnier, mensuel, ...) sur les données traitées
    - Variables utilisées par tout les modules
        -> zref
        -> filename
        -> ...
"""

def atlas_filename(diratlas, atlas_infos):
    """
    :param diratlas: Chemin vers le répertoire contenant les atlas
    :param atlas_infos: Dictionnaire contenant les champs suivants : (mode, date, reso, timeflag, typestat)
                      permettant de connaitre les informations sur les statistiques à calculer
    
    Fonction permettant de générer le nom d'un fichier atlas en fonction de ses
    paramètres contenus dans atlas_infos.
    
    :rtype: String
    """

    if False: # many subfolders - short name for the atlas
        atlas_name = '%s_%g_annual' % (atlas_infos['TYPESTAT'], atlas_infos['RESO'])
        ncfile = '%s/%s/%s/%s/%s/%s.nc' % (diratlas, atlas_infos['RESO'],
                                           atlas_infos['YEAR'], atlas_infos['MODE'],
                                           atlas_infos['TYPESTAT'], atlas_name)
        print(ncfile)

    else: # one folder - long name for the atlas
        ncfile ='%s/%s_%g_%s_%s.nc' % (diratlas, atlas_infos['TYPESTAT'], 
                                       atlas_infos['RESO'], atlas_infos['YEAR'], 
                                       atlas_infos['MODE'])
        print(ncfile)

    return ncfile

