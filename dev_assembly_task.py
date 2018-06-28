#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 08:08:58 2018

@author: therry
"""

import dev_zref as dev
import param as param
import tile as ti

import os
import pandas as pd

n_proc = 27
itile = 238

def synchronize_zref_profiles_from_task(itile, profiles, zref_profiles):
    """
    :param itile: Numéro de la tile que l'on souhaite générer
    :param argo_tile: Version antérieure de l'argo_%3i à synchroniser avec la nouvelle
    :param interp_result: Version mise à jour d'argo_tile avec compte rendu interpolation
    Fonction utilisée pour mettre à jour les argo_%3i avec les nouveaux profiles de la base Argo,
    mais également pour le retour d'information après l'interpolation des profiles

    :rtype: DataFrame
    """

    keys = ['CT', 'SA', 'RHO', 'BVF2']
    for tag in profiles['CT'].index:
        if tag in zref_profiles['CT'].index:
            for key in keys:
                print('Profile tagged %i changed' % tag)
                zref_profiles[key].loc[tag] = profiles[key].loc[tag]
    
    idx = profiles['CT'].index.difference(zref_profiles['CT'].index)
        
    for k in keys:
        zref_profiles[k] = zref_profiles[k].append(profiles[k].loc[idx])
    
    ti.write_zref_profiles(itile, zref_profiles)

    
keys = ['DATA_MODE', 'FLAG', 'JULD', 'LONGITUDE', 'LATITUDE', 'STATUS', 'WMO', 'IDAC', 'IPROF']
keys_prof = ['CT', 'SA', 'RHO', 'BVF2']
tile = pd.DataFrame(columns=keys)
CT = pd.DataFrame(columns = param.zref)
SA = pd.DataFrame(columns = param.zref)
RHO = pd.DataFrame(columns = param.zref)
BVF2 = pd.DataFrame(columns = param.zref)
profiles = {'CT' : CT, 'SA' : SA, 'RHO' : RHO, 'BVF2' :BVF2}

for i in range(n_proc):
    argo_task = dev.read_argo_task(itile, i)
    tile = pd.concat([tile, argo_task], sort=False)

    
for i in range(n_proc):
    task = dev.read_zref_task(itile, i)
    for key in keys_prof:
        profiles[key] = pd.concat([profiles[key], task[key]], sort=False)

argo_tile = ti.read_argo_tile(itile)
filename = '%s/zref_profiles_%003i.pkl' % (param.get_path('zref_profiles'), itile)
if os.path.isfile(filename):
    zref_profiles = ti.read_zref_profiles(itile)
else:
    CT = pd.DataFrame(columns = param.zref)
    SA = pd.DataFrame(columns = param.zref)
    RHO = pd.DataFrame(columns = param.zref)
    BVF2 = pd.DataFrame(columns = param.zref)
    zref_profiles = {'CT' : CT, 'SA' : SA, 'RHO' : RHO, 'BVF2' :BVF2}


ti.synchronize_argo_tile_from_interpolation(itile, argo_tile, tile)
zref_profiles = synchronize_zref_profiles_from_task(itile, profiles, zref_profiles)