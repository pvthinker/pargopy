#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 12:39:10 2018

@author: therry

Tools to generate and maintains the processed database

they are high-level routines that rely on smaller modules
"""

import pandas as pd

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
    argo_profile_dic = tools.read_profile(idac, wmo, data = True, header = True)
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
