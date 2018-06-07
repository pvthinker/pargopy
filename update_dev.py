#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 12:36:33 2018

@author: therry
"""
import pandas as pd
import numpy as np
import os
import random

import database as db
import general_tools as tools
import param as param

path_argo = "/home/ref-argo/gdac/dac"

daclist = ['aoml', 'bodc', 'coriolis', 'csio',
           'csiro', 'incois', 'jma', 'kma',
           'kordi', 'meds', 'nmdis']

def count_profiles_in_database(wmostats):
    """
    Count the total number of profiles in database
    :rtype: int
    """

    nbprofiles = 0
    for nbpr in wmostats['N_PROF']:
        nbprofiles += nbpr
    print('number of new profiles in Argo database: %i' % nbprofiles)
    return nbprofiles


def get_new_wmos_argo_global(wmos_infos):
    """Build argodb from the infos in wmos_infos

    Once it is created it is more efficient to read it from the disk
    using 'read_argodb()'
    
    :rtype: dic

    """
    keys = ['LONGITUDE', 'LATITUDE', 'JULD', 'DATA_MODE']
    keys_final = ['FLAG', 'LONGITUDE', 'LATITUDE', 'JULD', 'DATA_MODE', 'STATUS']
    key_int = ['TAG', 'FLAG', 'N_PROF']
    key_float = ['LONGITUDE', 'LATITUDE', 'JULD']
    key_char = ['DATA_MODE']
    key_bool = ['STATUS']

    n_profiles = count_profiles_in_database(wmos_infos)
    argodb = {}
    argo_wanted = {}

    for k in key_int:
        argodb[k] = np.zeros((n_profiles,), dtype=int)
    for k in key_float:
        argodb[k] = np.zeros((n_profiles,))
    for k in key_char:
        argodb[k] = np.zeros((n_profiles,), dtype='c')
    for k in key_bool:
        argodb[k] = np.zeros((n_profiles,), dtype=bool)

    iprof = 0

    for idx in wmos_infos.index:
        kdac = wmos_infos.loc[idx, 'DACID']
        wmo = idx # wmos are index of wmos_infos
        dac = param.daclist[kdac]
        output = db.read_profile(dac, wmo, header=True, verbose=False)
        n_prof = output['N_PROF']

        for key in keys:
            argodb[key][iprof:iprof+n_prof] = output[key]
        for kprof in range(n_prof):
            tag = tools.get_tag(kdac, wmo, kprof)
            argodb['TAG'][iprof+kprof] = tag

        iprof += n_prof
        print(' {:8,}/{:,} : {} - {}'.format(iprof, n_profiles, dac, wmo))

    for k in keys_final:
        argo_wanted[k] = argodb[k]

    argo_new_wmos = pd.DataFrame(argo_wanted, index = argodb['TAG'])

    return argo_new_wmos

# On récupère l'argo_global existant ( vide ou non)
argo_global = db.read_argo_global()
#==============================================================================
# argo_global['STATUS'] = True
# list_flag = random.sample(range(1800000), 200000)
# argo_global.loc[argo_global.index[list_flag], 'FLAG'] = 111
#==============================================================================
wmodic = db.create_wmodic()

# On récupère l'ancien wmos_infos si il existe, sinon on le crée vide
if os.path.exists('%s/wmos_infos_2017.pkl' % param.get_path('database')):
    old_wmos_infos = pd.read_pickle('%s/wmos_infos_2017.pkl' % param.get_path('database'))
else:
    old_wmos_infos = pd.DataFrame(columns = ['DACID', 'DATE_UPDATE'])

# On crée le nouveau wmos_infos
#  db.get_wmo_infos(param.get_path('argo'))
wmos_infos = db.read_wmos_infos()

##############################################################################
#                       Récupération des nouveaux wmos
##############################################################################

#  On récupère les index des nouveaux wmos (WMO)
idx = wmos_infos.index.difference(old_wmos_infos.index)

# Avec ces index, on récupère le TAG correspondant
counter = 0
for wmo in idx:
    for kprof in range(wmos_infos.loc[wmo, 'N_PROF']):
        counter += 1

new_tags_list = [0 for i in range(counter)]

i = 0

# On crée une liste contenant les nouveaux tags
for wmo in idx:
    for kprof in range(wmos_infos.loc[wmo, 'N_PROF']):
        tag = tools.get_tag(wmos_infos.loc[wmo, 'DACID'], wmo, kprof)
        new_tags_list[i] = tag
        i += 1

# On crée un sous ensemble d'argo_global contenant les nouveaux wmos
new_wmos_argo = get_new_wmos_argo_global(wmos_infos.loc[idx])


##############################################################################
#                       Mise à jour des wmos modifiés
##############################################################################
modified_tag_list = []
new_prof_tags_list = []
nb_new_prof_list = []
idx_new_prof = []
idx_unchanged_prof = []
for wmo in wmos_infos.index:
    if wmo in old_wmos_infos.index:
        # si le wmo est dans old, on regarde si la valeur d'une des valeurs de la ligne associée à change
        changed = any(wmos_infos.loc[wmo, :] != old_wmos_infos.loc[wmo, :])
        if changed:
            # Si les données ont changées, on regarde si le nombre de profile à augmenté
            if wmos_infos.loc[wmo, 'N_PROF'] == old_wmos_infos.loc[wmo, 'N_PROF']:
                # Si le nombre de profile n'a pas changé, on récupère les lignes
                # associées dans argo_global et on change le champ STATUS à False
                for kprof in range(wmos_infos.loc[wmo, 'N_PROF']):
                    tag = tools.get_tag(wmos_infos.loc[wmo, 'DACID'], wmo, kprof)
                    modified_tag_list.append(tag)
                
            else:
                # Si le nombre de profile a augmenté, il faut ajouter les nouveaux
                # TAG dans argo_global
                # on récupère le nombre de nouveaux profiles de chaque wmo
                nb_new_prof_list.append(wmos_infos.loc[wmo, 'N_PROF'] - old_wmos_infos.loc[wmo, 'N_PROF'])
                # on récupère les wmos
                idx_new_prof.append(wmo)

        else:
            # si aucune valeur correspondant aux wmos n'a changé,
            # On recopie simplement les lignes associées dans un sous ensemble
            # d'argo_global
            for kprof in range(wmos_infos.loc[wmo, 'N_PROF']):
                    tag = tools.get_tag(wmos_infos.loc[wmo, 'DACID'], wmo, kprof)
                    idx_unchanged_prof.append(tag)
    else:
        # si le wmo n'est pas dans old_wmos_infos
        # on crée un sous ensemble d'argo_global avec ses profiles
        pass

# On crée un autre sous ensemble d'argo_global pour les profiles modifiés
modified_prof_argo = argo_global.loc[modified_tag_list]
modified_prof_argo.loc[:, 'STATUS'] = False


# On crée un deuxième sous-ensemble argo_global pour les nouveaux profiles ajoutés
# et ceux qui ont été modifiés au sein de ces wmo
new_prof_argo = get_new_wmos_argo_global(wmos_infos.loc[idx_new_prof])

# on récupère les nouveaux TAG apparus dans argo_extraction
new_idx = new_prof_argo.index.difference(argo_global.index)

# on récupère les autres indices d'argo_extraction pour le traitement de ces données
idx_old_prof = [k for k in new_prof_argo.index if k not in new_idx]
new_prof_idx = []
for i in idx_old_prof:
    if new_prof_argo.loc[i, 'DATA_MODE'] == argo_global.loc[i, 'DATA_MODE']:
        new_prof_idx.append(i)
    else:
        pass

new_prof_argo.loc[new_prof_idx] = argo_global.loc[new_prof_idx]

# on recopie les autres profiles au sein d'argo_global comme ils n'ont pas changé
unchanged_prof_argo = argo_global.loc[idx_unchanged_prof]

##############################################################################
#                       Assemblage d'argo_global mis à jour
##############################################################################

argo_global_update = pd.concat([new_wmos_argo, modified_prof_argo, new_prof_argo, unchanged_prof_argo])
# Une fois tout les sous ensembles assemblés, on écrit la nouvelle version
# d'argo_global dans son fichier pickle
db.write_argo_global(argo_global_update)
argo_global = pd.read_pickle('%s/argo_global.pkl' % param.get_path('database'))
argo_global_2017 = pd.read_pickle('%s/argo_global_2017.pkl' % param.get_path('database'))
print('-'*60)
print('Quantity of modified tag : %i' % len(modified_tag_list))
print('Quantity of unchanged profiles in wmos with new prof : %i' % len(new_prof_idx))
print('Quantity of unchanged profiles in unchanged wmos : %i' % len(idx_unchanged_prof))
y = 0
for i in nb_new_prof_list:
    y += i
print('Quantity of new_prof : %i' % (y+len(new_tags_list)))
print('-'*60)