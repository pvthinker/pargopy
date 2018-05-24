#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 06:52:59 2018

@author: therry

Ideas for function melting
"""

import pickle as pickle
#  import param as param
#  path_localdata = param.path_to_data

#  ----------------------------------------------------------------------------
def write_dic(name, dic, path_localdata):
    """
    Function used to write each dic used to create .pkl files
    Regroups :
    - write_wmodic, write_wmstats, write_argodb from argodb.py
    - write_argo_filter from research_tools.py
    - write_tile from tile.py
    
    :rtype: None
    """

    with open('%s/%s.pkl' % (path_localdata, name), 'w') as f:
        pickle.dump(dic, f)

#  ----------------------------------------------------------------------------
def read_dic(name, path_localdata):
    """
    Function used to read each dic used to create .pkl files
    Regroups :
    - read_wmodic, read_wmstats, read_argodb from argodb.py
    - read_argo_filter from research_tools.py
    - read_tile from tile.py
    
    :rtype: dict"""

    print('read %s.pkl' % name)
    with open('%s/%s.pkl' % (path_localdata, name), 'r') as f:
        dic = pickle.load(f)
    return dic