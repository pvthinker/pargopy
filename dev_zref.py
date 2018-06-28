#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 07:26:36 2018

@author: therry
"""

import pandas as pd

import param as param
import tile as ti
import interpolation_tools as interp
import general_tools as tools



def read_argo_task(itile, itask):
    """
    Fonction utilisée pour récupérer la DataFrame 'argo_global' du fichier pickle
    dans lequel elle a été préalablement sauvegardé.

    :rtype: DataFrame
    """
    argo_task = pd.read_pickle('%s/task_%003i_%003i.pkl' % (param.get_path('task'), itile, itask))

    return argo_task


def write_argo_task(itile, argo_task, itask):
    """
    :param itile: Numéro de la tile que l'on veut sauvegarder
    :param argo_tile: Tile à sauvegarder

    Fonction utilisée pour écrire notre DataFrame argo_tile au sein d'un fichier
    pickle portant le nom de : argo_%003i.pkl.

    :rtype: None
    """
    argo_task.to_pickle('%s/task_%003i_%003i.pkl' % (param.get_path('task'), itile, itask))


def write_zref_task(itile, zref_task, itask):
    """
    :param itile: Numéro de la tile que l'on veut sauvegarder
    :param zref_profiles: Dictionnaire contenant les profiles interpolés

    Fonction utilisée pour écrire nos profiles interpolés dans un fichier pickle

    :rtype: None
    """
    filename = '%s/task_%003i_%003i.pkl' % (param.get_path('zref_task'), itile, itask)
    pd.to_pickle(zref_task, filename)


def read_zref_task(itile, itask):
    """
    :param itile: Numéro de la tile que l'on veut sauvegarder
    :param zref_profiles: Dictionnaire contenant les profiles interpolés

    Fonction utilisée pour écrire nos profiles interpolés dans un fichier pickle

    :rtype: None
    """
    filename = '%s/task_%003i_%003i.pkl' % (param.get_path('zref_task'), itile, itask)
    print('Read zref task: %s' % filename)
    dic = pd.read_pickle(filename)

    return dic


def generate_argo_task(tile, itile):
    """
    Fonction divisant le travail a effectuer en fonction du nombre de processeur
    """
    
    #  keys = ['DATA_MODE', 'FLAG', 'JULD', 'LONGITUDE', 'LATITUDE', 'STATUS', 'WMO', 'IDAC', 'IPROF']
    #  argodb = pd.DataFrame(columns=keys)
    
    n_proc = 27
    
    subargo = tile[tile['STATUS'] == False]
    all_tags_infos = tools.retrieve_infos_from_tag(subargo.index)
    subargo.loc[:,'WMO'] = all_tags_infos['WMO']
    subargo.loc[:,'IDAC'] = all_tags_infos['IDAC']
    subargo.loc[:,'IPROF'] = all_tags_infos['IPROF']
    
    n_task = len(subargo.index) / n_proc
    
    i0 = 0
    i1 = 0
    
    for i in range(n_proc):
        i1 = i0 + n_task
        if i == (n_proc - 1):
            argo_task = subargo.iloc[i0:]
        else:
           argo_task = subargo.iloc[i0:i1]
        write_argo_task(itile, argo_task, i) 
        i0 += n_task
    
    
#==============================================================================
#     for i in range(n_proc):
#         argo_task = read_argo_task(itile, i)
#         argodb = pd.concat([argodb, argo_task], sort=False)
#==============================================================================


def generate_zref_task(task, itask, itile):
    """
    Fonction divisant le travail a effectuer en fonction du nombre de processeur
    """
    
    keys = ['CT', 'SA', 'RHO', 'BVF2']
    
    if len(task.index) == 0:
        zref_prof_to_update = []
    else:
        zref_profiles, interp_result = interp.interpolate_profiles(task)
    
    CT = pd.DataFrame(columns = param.zref)
    SA = pd.DataFrame(columns = param.zref)
    RHO = pd.DataFrame(columns = param.zref)
    BVF2 = pd.DataFrame(columns = param.zref)
    zref_prof_to_update = {'CT' : CT, 'SA' : SA, 'RHO' : RHO, 'BVF2' :BVF2}
    
    for k in keys:
            zref_prof_to_update[k] = zref_prof_to_update[k].append(zref_profiles[k])
    write_zref_task(itile, zref_prof_to_update, itask)
    print('zref_task_%003i_%003i generated' % (itile, itask))


def main(itile, itask):
    tile = ti.read_argo_tile(itile)
    # on peut lui envoyer une tile ou une portion d'argo_global, ou encore plusieurs
    # tiles concatenees au préalable ou encore argo_global entier
    generate_argo_task(tile, itile)
    task = read_argo_task(itile, itask)
    generate_zref_task(task, itask, itile)

if __name__ == '__main__':

    n_proc = 27
    itile = 238
    main(itile, 1)
