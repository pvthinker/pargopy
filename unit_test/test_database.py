#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:53:00 2018

@author: therry

Module contenant la batterie de tests unitaires pour vérifier le fonctionnement
des fonctions de general_tools.py

"""
import os
import pytest
import sys
import pandas as pd
sys.path[:0] = ['../']
import database as db
import param

@pytest.fixture()
def parameters_definer():
    idac = 0
    wmo = 1900060
    iprof= 0
    argo_profile_dic = db.read_profile(idac, wmo, data = True, header = True)
    return (idac, wmo, iprof, argo_profile_dic)


def test_get_prof(parameters_definer):
    idac, wmo, iprof, argo_profile_dic = parameters_definer
    data = db.get_prof(iprof, argo_profile_dic)
    assert 'TEMP' in data.keys()
    assert 'PSAL' in data.keys()


def test_argo_profile_dic_to_dataframe(parameters_definer):
    idac, wmo, iprof, argo_profile_dic = parameters_definer
    argo_profile_dataframe = db.argo_profile_dic_to_dataframe(idac, wmo)
    keys = ['DATA_MODE', 'LONGITUDE', 'LATITUDE', 'JULD', 'DATE_UPDATE', 'FLAG', 'PROFILE_DATA']
    for k in keys:
        assert k in argo_profile_dataframe.keys()

def test_create_workspace():
    paths_list = ['pargopy_output', 'database', 'parallel', 'zref_profiles', 'stats', 'atlas']
    for path in paths_list:
        assert os.path.exists(param.get_path(path))
    argo_global = pd.read_pickle('%s/argo_global.pkl' % param.get_path('database'))
    keys = ['DATA_MODE', 'LONGITUDE', 'LATITUDE', 'JULD', 'DATE_UPDATE', 'FLAG', 'PROFILE_DATA']
    for k in keys:
        assert k in argo_global.keys()