#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:53:00 2018

@author: therry

Module contenant la batterie de tests unitaires pour vérifier le fonctionnement
des fonctions de general_tools.py

"""
import pytest
import sys
sys.path[:0] = ['../']
import database as db
import general_tools as tools

@pytest.fixture()
def parameters_definer():
    idac = 0
    wmo = 1900060
    iprof= 0
    argo_profile_dic = tools.read_profile(idac, wmo, data = True, header = True)
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