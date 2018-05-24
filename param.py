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

