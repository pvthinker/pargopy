#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 12:39:10 2018

@author: therry

Tools to generate and maintains the processed database

they are high-level routines that rely on smaller modules
"""

import os
import pandas as pd
from netCDF4 import Dataset
import numpy as np

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
                 verbose=True):
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

    argo_profile_dic = {}

    required_keys = set(['TEMP', 'PSAL', 'PRES'])

    if (os.path.isfile(filename)):
        with Dataset(filename, "r", format="NETCDF4") as f:
            argo_profile_dic['DACID'] = param.daclist.index(dac)
            argo_profile_dic['WMO'] = wmo
            argo_profile_dic['N_PROF'] = len(f.dimensions['N_PROF'])
            argo_profile_dic['N_LEVELS'] = len(f.dimensions['N_LEVELS'])
            # DATE_UPDATE is an array of 14 characters in the *_prof.nc
            # we transform it into an int
            # YYYYMMDDhhmmss
            argo_profile_dic['DATE_UPDATE'] = ''.join(f.variables['DATE_UPDATE'][:])

            keyvar = set(f.variables.keys())

            if required_keys.issubset(keyvar):
                argo_profile_dic['TSP_QC'] = '1'
            else:
                argo_profile_dic['TSP_QC'] = '2'

            if header or headerqc or data or dataqc:
                if iprof is None:
                    idx = range(argo_profile_dic['N_PROF'])
                    argo_profile_dic['IPROF'] = np.arange(argo_profile_dic['N_PROF'])
                else:
                    idx = iprof
                    argo_profile_dic['IPROF'] = iprof

            if header:
                for key in key_header:
                    argo_profile_dic[key] = f.variables[key][idx]
                    argo_profile_dic['DATA_MODE'] = np.asarray(
                        [c for c in f.variables['DATA_MODE'][idx]])

            if headerqc:
                for key in key_headerqc:
                    argo_profile_dic[key] = f.variables[key][idx]

            if data:
                for key in key_data:
                    if argo_profile_dic['TSP_QC'] == '1':
                        argo_profile_dic[key] = f.variables[key][idx, :]
                    else:
                        argo_profile_dic[key] = np.NaN+np.zeros(
                            (argo_profile_dic['N_PROF'], argo_profile_dic['N_LEVELS']))

            if dataqc:
                for key in key_dataqc:
                    if argo_profile_dic['TSP_QC'] == '1':
                        argo_profile_dic[key] = f.variables[key][idx]
                    else:
                        argo_profile_dic[key] = np.zeros(
                            (argo_profile_dic['N_PROF'], argo_profile_dic['N_LEVELS']), dtype=str)

    return argo_profile_dic


def argo_profile_dic_to_dataframe(idac, wmo):
    """
    :param kdac: Index of the dac (aoml = 1, bodc = 2, coriolis = 3, ...)
    :param wmo: WMO number
    Fonction retournant les profiles argo sous forme de DataFrame.
    Le TAG associé au profile est passé en index des colonnes suivantes :
        - DATA_MODE
        - LONGITUDE
        - LATITUDE
        - JULD
        - DATE_UPDATE
        - PROFILE_DATA
        - FLAG
    
    La pression (PRES) est passée en index des colonnes suivantes :
        - TEMP
        - PSAL
    
    TEMP et PSAL font partie d'un DataFrame contenu dans PROFILE_DATA.

    :rtype: DataFrame
    """
    argo_profile_dic = read_profile(idac, wmo, data = True, header = True)
    keys = ['DATA_MODE', 'LONGITUDE', 'LATITUDE', 'JULD', 'DATE_UPDATE']
    
    argo_profile_dic_V2 = {}
    # assign a tag to each profile
    tag = [tools.get_tag(idac, wmo, iprof) for iprof in range(argo_profile_dic['N_PROF'])]
    # set flag to 0 by default, 0 = good profile
    argo_profile_dic_V2['FLAG'] = [0 for iprof in range(argo_profile_dic['N_PROF'])]
    # the date of last update is the same for all the profiles of a wmo
    argo_profile_dic_V2['DATE_UPDATE'] = [argo_profile_dic['DATE_UPDATE'] for iprof in range(argo_profile_dic['N_PROF'])]

    for k in keys:
        argo_profile_dic_V2[k] = argo_profile_dic[k]

    # si on veut aussi retourner les 'data' = les profils
    # 1 profile (T,S,p) = 1 DataFrame
    argo_profile_dic_V2['PROFILE_DATA'] = [get_prof(iprof, argo_profile_dic) for iprof in range(argo_profile_dic['N_PROF'])]

    argo_profile_dataframe = pd.DataFrame(argo_profile_dic_V2, index = tag)

    return argo_profile_dataframe


def get_prof(iprof, argo_profile_dic):
    """
    :param iprof: Index of the profile
    :param argo_profile_dic: Dictionnary containing the Argo profiles

    Transform profile 'iprof' stored in 'res' into a DataFrame
    whose index is 'PRES' and the two columns are 'TEMP' and 'PSAL'
    also removes all NaN
    
    :rtype: DataFrame
    """
    mskT = argo_profile_dic['TEMP'][iprof, :].mask 
    mskS = argo_profile_dic['PSAL'][iprof, :].mask 
    mskP = argo_profile_dic['PRES'][iprof, :].mask 
    # mskT == True means the data should be omitted
    # msk = True means T, S and P are ok
    msk = ~(mskT | mskS | mskP)

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
        - DATE_UPDATE
        - PROFILE_DATA
        - FLAG
    Ce DataFrame sera ensuite sauvegardé au sein d'un fichier pickle nommé argodb.pkl
    
    :rtype: None
    """
    keys = ['DATA_MODE', 'LONGITUDE', 'LATITUDE', 'JULD', 'DATE_UPDATE', 'FLAG', 'PROFILE_DATA']
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


def update_argo_global():
    """
    Fonction utilisée pour mettre à jour argo_global
    """



