#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 12:39:10 2018

@author: therry

Tools to generate and maintains the processed database

they are high-level routines that rely on smaller modules
"""

import os
import glob
import pandas as pd
from netCDF4 import Dataset
import numpy as np
import time as time

import param as param
import general_tools as tools

header_keys = ['DATA_MODE', 'LONGITUDE', 'LATITUDE', 'JULD', 'FLAG', 'STATUS']


def synchronize_headers():
    """ 
    Propagate the header of each tile onto the global header

    """
    header_global = tools.create_empty_header()
    # pour crÃ©er un empty DataFrame avec les bonnes clefs:
    header_global = pd.DataFrame(columns=header_keys+['TAG'])
    header_global = header_global.set_index('TAG')
    
    for itile in range(ntiles):
        header_tile = tile.read_header(itile)
        header_global =header_global.merge(header_tile,
                                           left_index=True, right_index=True,
                                           how='outer', on=header_keys)
        # warning: because of the overlapping of tiles, some tags
        # appear in several tiles, they need to be exactly identical
        # otherwise merge will complain. In that case, it means that
        # we do someting wrong
    tools.write_header_global(header_global)

def read_profile(dac, wmo, iprof=None,
                 header=False, data=False,
                 headerqc=False, dataqc=False,
                 shortheader=False,
                 verbose=True, path=None):
    """
    :param dac: DAC du profil recherché
    :param wmo: WMO du profil recherché
    :param iprof: Numéro du profil recherché
    :param header: Sélectionne seulement LATITUDE, LONGITUDE et JULD
    :param headerqc: Sélectionne seulement POSITION_QC et JULD_QC
    :param data: Sélectionne TEMP, PSAL et PRES
    :param dataqc: Sélectionne TEMP_QC, PSAL_QC et PRES_QC
    :param verbose: ???
    
    Les valeurs sélectionnée grâce aux arguments passés à la fonction définissent
    la DataFrame que retournera celle-ci.
    
    Basic driver to read the \*_prof.nc data file

    The output is a dictionnary of vectors
    - read one or all profiles read the header (lat, lon, juld) or not
    - read the data or not always return IDAC, WMO, N_PROF, N_LEVELS
    - and DATA_UPDATE (all 5 are int)

    :rtype: dict
    """
    key_header = ['LATITUDE', 'LONGITUDE', 'JULD']
    key_headerqc = ['POSITION_QC', 'JULD_QC']
    key_data = ['TEMP', 'PSAL', 'PRES']
    key_dataqc = ['TEMP_QC', 'PSAL_QC', 'PRES_QC']

    if type(dac) is int:
        dac = param.daclist[dac]

    filename = param.get_argo_filename(dac, wmo)

    if verbose:
        print('/'.join(filename.split('/')[-3:]))

    output = {}

    required_keys = set(['TEMP', 'PSAL', 'PRES'])

    if (os.path.isfile(filename)):
        with Dataset(filename, "r", format="NETCDF4") as f:
            output['DACID'] = param.daclist.index(dac)
            output['WMO'] = wmo
            output['N_PROF'] = len(f.dimensions['N_PROF'])
            output['N_LEVELS'] = len(f.dimensions['N_LEVELS'])
            # DATE_UPDATE is an array of 14 characters in the *_prof.nc
            # we transform it into an int
            # YYYYMMDDhhmmss
            output['DATE_UPDATE'] = ''.join(f.variables['DATE_UPDATE'][:])
            if shortheader:
                pass
            else:
                keyvar = set(f.variables.keys())

                if required_keys.issubset(keyvar):
                    output['TSP_QC'] = '1'
                else:
                    output['TSP_QC'] = '2'

                if header or headerqc or data or dataqc:
                    if iprof is None:
                        idx = range(output['N_PROF'])
                        output['IPROF'] = np.arange(output['N_PROF'])
                    else:
                        idx = iprof
                        output['IPROF'] = iprof

                if header:
                    for key in key_header:
                        output[key] = f.variables[key][idx]
                        output['DATA_MODE'] = np.asarray(
                            [c for c in f.variables['DATA_MODE'][idx]])

                if headerqc:
                    for key in key_headerqc:
                        output[key] = f.variables[key][idx]

                if data:
                    for key in key_data:
                        if output['TSP_QC'] == '1':
                            output[key] = f.variables[key][idx, :]
                        else:
                            output[key] = np.NaN+np.zeros(
                                (output['N_PROF'], output['N_LEVELS']))

                if dataqc:
                    for key in key_dataqc:
                        if output['TSP_QC'] == '1':
                            output[key] = f.variables[key][idx]
                        else:
                            output[key] = np.zeros((output['N_PROF'],
                                                    output['N_LEVELS']),
                                                   dtype=str)
    return output


def get_new_wmos_argo_global(wmos_infos):
    """
    :param wmos_infos: Dataframe contenant les informations d'un wmo
    Fonction retournant les profiles argo sous forme de DataFrame.
    Le TAG associé au profile est passé en index des colonnes suivantes :
        - DATA_MODE
        - LONGITUDE
        - LATITUDE
        - JULD
        - DATE_UPDATE
        - FLAG

    :rtype: DataFrame
    """
    keys = ['LONGITUDE', 'LATITUDE', 'JULD', 'DATA_MODE']
    keys_final = ['FLAG', 'LONGITUDE', 'LATITUDE', 'JULD', 'DATA_MODE', 'STATUS']
    key_int = ['TAG', 'FLAG', 'N_PROF']
    key_float = ['LONGITUDE', 'LATITUDE', 'JULD']
    key_char = ['DATA_MODE']
    key_bool = ['STATUS']

    n_profiles = tools.count_profiles_in_database(wmos_infos)
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
        output = read_profile(dac, wmo, header=True, verbose=False)
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


def get_prof(iprof, argo_profile_dic):
    """
    :param iprof: Index of the profile
    :param argo_profile_dic: Dictionnary containing the Argo profiles

    Transform profile 'iprof' stored in 'res' into a DataFrame
    whose index is 'PRES' and the two columns are 'TEMP' and 'PSAL'
    also removes all NaN
    
    :rtype: DataFrame
    """
    if (type(argo_profile_dic['TEMP'][iprof, :]) == np.ma.core.MaskedArray and 
        type(argo_profile_dic['PSAL'][iprof, :]) == np.ma.core.MaskedArray and 
        type(argo_profile_dic['PRES'][iprof, :]) == np.ma.core.MaskedArray):
        mskT = argo_profile_dic['TEMP'][iprof, :].mask 
        mskS = argo_profile_dic['PSAL'][iprof, :].mask 
        mskP = argo_profile_dic['PRES'][iprof, :].mask 
    # mskT == True means the data should be omitted
    # msk = True means T, S and P are ok
    else:
        mskT = False
        mskS = False 
        mskP = False

    msk = ~(mskT | mskS | mskP)

    if type(msk) != np.ndarray:
        pass
    elif len(msk) <= 5:
        pass
    else:
        profile_data = pd.DataFrame({'TEMP': argo_profile_dic['TEMP'][iprof][msk],
                                     'PSAL': argo_profile_dic['PSAL'][iprof][msk]},
                                    index=argo_profile_dic['PRES'][iprof][msk])

        return profile_data


def create_database():
    """
    Fonction à utiliser si argodb.pkl n'existe pas (lors d'une première utilisation
    du programme par exemple).
    Cette fonction crée une DataFrame vide contenant simplement les clefs des futures
    variables qui y seront ajoutées :
        - DATA_MODE
        - LONGITUDE
        - LATITUDE
        - JULD
        - STATUS
        - FLAG
    Ce DataFrame sera ensuite sauvegardé au sein d'un fichier pickle nommé argodb.pkl
    
    :rtype: None
    """
    keys = ['DATA_MODE', 'LONGITUDE', 'LATITUDE', 'JULD', 'STATUS', 'FLAG']
    argodb = pd.DataFrame(columns = keys)
    argodb.to_pickle('%s/argo_global.pkl' % param.get_path('database'))
    for i in range(300):
            argodb.to_pickle('%s/argo_%003i.pkl' % (param.get_path('parallel'), i))


def create_workspace():
    """
    Fonction utilisée pour mettre en place l'espace de travail lors de la 
    première utilisation du programme :
        - Création des différents répertoires utilisés pour le stockage des fichiers
          générés par le programme lors de son execution
        - Création de nos fichiers pickle stockés dans le répertoire 'database' situé
          dans notre workdir.
          create_argodb génère un fichier destiné à contenir un snapchot d'Argo global, ainsi
          que 300 fichiers pickle destinés à diviser cette base de donnée en 300 dalles
          utilisées ultérieurement pour le calcul parallèle
    
    :rtype: None
    """
    create_workdir()
    create_database()

def create_workdir():
    """
    Fonction permettant de créer les différents répertoires utilisés pour l'espace
    de travail :
        - pargopy_output/
            - database/
                - parallel/
            - zref_profiles/
            - stats_work/
            - atlas/
                - reso_1.0/
                - reso_0.5/
                - reso_0.25/
    
    :rtype: None
    """
    paths_list = ['pargopy_output', 'database', 'parallel', 'zref_profiles', 'stats', 'atlas', 'reso']
    reso_list = [0.25, 0.5, 1.0]

    for path in paths_list:
        if path == 'reso':
            for reso in reso_list:
                os.mkdir(param.get_path(path, reso))
                print('File created : %s' % param.get_path(path, reso))
        else:
            os.mkdir(param.get_path(path))
            print('File created : %s' % param.get_path(path))


def write_argo_global(argo_global):
    """
    :param argo_global: DataFrame à écrire au sein du fichier pickle

    Fonction utilisée pour écrire notre DataFrame argo_global au sein d'un fichier
    pickle portant le nom de : argo_global.pkl.

    :rtype: None
    """
    argo_global.to_pickle('%s/argo_global.pkl' % param.get_path('database'))


def read_argo_global():
    """
    Fonction utilisée pour récupérer la DataFrame 'argo_global' du fichier pickle
    dans lequel elle a été préalablement sauvegardé.

    :rtype: DataFrame
    """
    argo_global = pd.read_pickle('%s/argo_global.pkl' % param.get_path('database'))

    return argo_global


def get_wmo_infos(path_gdac):
    w = []
    d = []
    n = []
    u = []
    for dac in param.daclist:  # ['coriolis']
        prfiles = glob.glob('{}/{}/*/*_prof.nc'.format(path_gdac, dac))
        wmos = [int(f.split('/')[-2]) for f in prfiles]
        for wmo in wmos:
            output = read_profile(dac, wmo, shortheader=True, verbose=True,
                                     path=path_gdac)
            w += [wmo]
            d += [param.daclist.index(dac)]
            n += [output['N_PROF']]
            u += [output['DATE_UPDATE']]

    wmos_infos = pd.DataFrame({'DACID': d, 'N_PROF': n, 'DATE_UPDATE': u}, index=w)

    wmos_infos.to_pickle('%s/wmos_infos.pkl' % param.get_path('database'))


def flag_argodb(argo_global, wmodic):
    """
    Add the flag to argodb

    :rtype: dic
    """

    tag = argo_global.index[argo_global['FLAG'] == 0]
        

    infos = tools.retrieve_infos_from_tag(tag)
    wmos = set(infos['WMO'])
    idx = []
    tag_101 = []
    tag_102 = []
    tag_103 = []
    for w in wmos:
        print('add flag to wmo %i' % w)
        idx = np.where(infos['WMO'] == w)[0]
        iprof = infos['IPROF'][idx]
        dac = tools.dac_from_wmo(wmodic, w)
        res = read_profile(dac, w, headerqc=True)
        # print('-'*80)
        # print(res['TEMP_QC'].data)
        # print(res['POSITION_QC'].data)
        for k in iprof:
            if res['POSITION_QC'][k] != '1':
                 tag_101.append(tools.get_tag(param.daclist.index(dac), w, k))
            if res['JULD_QC'][k] != '1':
                tag_102.append(tools.get_tag(param.daclist.index(dac), w, k))
            if res['TSP_QC'] != '1':
                tag_103.append(tools.get_tag(param.daclist.index(dac), w, k))
    
    argo_global.loc[tag_101, 'FLAG'] = 101
    argo_global.loc[tag_102, 'FLAG'] = 102
    argo_global.loc[tag_103, 'FLAG'] = 103
                   
    return argo_global


def read_wmos_infos():
    """
    Fonction utilisée pour récupérer la DataFrame 'wmos_infos' du fichier pickle
    dans lequel elle a été préalablement sauvegardé.

    :rtype: DataFrame
    """
    wmos_infos = pd.read_pickle('%s/wmos_infos.pkl' % param.get_path('database'))

    return wmos_infos


#  ----------------------------------------------------------------------------
def create_wmodic():
    """
    Fonction permettant de créer un dictionnaire contenant la liste des dac et
    des wmos
    
    :rtype: dict
    """
    wmodic = {}
    for idac, dac in enumerate(param.daclist):
        prfiles = glob.glob('{}/{}/*/*_prof.nc'.format(param.get_path('argo'), dac))
        wmodic[dac] = [int(f.split('/')[-2]) for f in prfiles]
       
    return wmodic


def update_argo_global():
    """
    Fonction utilisée pour mettre à jour argo_global avec la base de donnée Argo
    évolutive.
    La mise à jour a lieu si :
        - l'attribut DATE_UPDATE contient une date plus récente que celle contenue
          dans argo_global
        - le TAG obtenu en manipulant le dac, wmo, iprof d'un profile est inconnu
          dans argo_global

    :rtype: DataFrame
    """
    wmodic = create_wmodic()
    argo_global = read_argo_global()

    if os.path.exists('%s/wmos_infos.pkl' % param.get_path('database')):
        old_wmos_infos = pd.read_pickle('%s/wmos_infos_2017.pkl' % param.get_path('database'))
    else:
        old_wmos_infos = pd.DataFrame(columns = ['DACID', 'DATE_UPDATE'])
    
    get_wmo_infos(param.get_path('argo'))
    wmos_infos = read_wmos_infos()

    # Récupération des nouveaux wmos

    idx = wmos_infos.index.difference(old_wmos_infos.index)

    counter = 0
    for wmo in idx:
        for kprof in range(wmos_infos.loc[wmo, 'N_PROF']):
            counter += 1
            
            new_tags_list = [0 for i in range(counter)]
    for wmo in idx:
        for kprof in range(wmos_infos.loc[wmo, 'N_PROF']):
            tag = tools.get_tag(wmos_infos.loc[wmo, 'DACID'], wmo, kprof)
            new_tags_list[i] = tag

    new_wmos_argo = get_new_wmos_argo_global(wmos_infos.loc[idx])


    modified_tag_list = []
    nb_new_prof_list = []
    idx_new_prof = []
    idx_unchanged_prof = []
    for wmo in wmos_infos.index:
        if wmo in old_wmos_infos.index:
            changed = any(wmos_infos.loc[wmo, :] != old_wmos_infos.loc[wmo, :])
            if changed:
                if wmos_infos.loc[wmo, 'N_PROF'] == old_wmos_infos.loc[wmo, 'N_PROF']:
                    for kprof in range(wmos_infos.loc[wmo, 'N_PROF']):
                        tag = tools.get_tag(wmos_infos.loc[wmo, 'DACID'], wmo, kprof)
                        modified_tag_list.append(tag)
                
                else:
                    nb_new_prof_list.append(wmos_infos.loc[wmo, 'N_PROF'] - old_wmos_infos.loc[wmo, 'N_PROF'])
                    idx_new_prof.append(wmo)

            else:
                for kprof in range(wmos_infos.loc[wmo, 'N_PROF']):
                    tag = tools.get_tag(wmos_infos.loc[wmo, 'DACID'], wmo, kprof)
                    idx_unchanged_prof.append(tag)
        else:
            pass

    modified_prof_argo = argo_global.loc[modified_tag_list]
    modified_prof_argo.loc[:, 'STATUS'] = False

    new_prof_argo = get_new_wmos_argo_global(wmos_infos.loc[idx_new_prof])

    new_idx = new_prof_argo.index.difference(argo_global.index)

    idx_old_prof = [k for k in new_prof_argo.index if k not in new_idx]
    new_prof_idx = []
    for i in idx_old_prof:
        if new_prof_argo.loc[i, 'DATA_MODE'] == argo_global.loc[i, 'DATA_MODE']:
            new_prof_idx.append(i)
        else:
            pass

    new_prof_argo.loc[new_prof_idx] = argo_global.loc[new_prof_idx]

    unchanged_prof_argo = argo_global.loc[idx_unchanged_prof]

    # Assemblage des sous ensembles
    argo_global_update = pd.concat([new_wmos_argo, modified_prof_argo, new_prof_argo, unchanged_prof_argo])
    argo_global_update = flag_argodb(argo_global_update, wmodic)
    write_argo_global(argo_global_update)
    print('Quantity of modified tag : %i' % len(modified_tag_list))
    print('Quantity of unchanged profiles in wmos with new prof : %i' % len(new_prof_idx))
    print('Quantity of unchanged profiles in unchanged wmos : %i' % len(idx_unchanged_prof))
    z = 0
    for i in nb_new_prof_list:
        z += i
    print('Quantity of new_prof : %i' % (z+len(new_tags_list)))