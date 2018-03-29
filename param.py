# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:19:12 2018

@author: herry
"""


#  ----------------------------------------------------------------------------
def data_argodb():
    path_to_argolog = '/local/tmp/1/herry/data/argo_log.nc'
    return path_to_argolog


#  ----------------------------------------------------------------------------
def data_research_tools():
    path_to_argotiles = '/local/tmp/1/herry/filter/argo_tile'
    path_to_argolog = data_argodb
    return path_to_argotiles, path_to_argolog


#  ----------------------------------------------------------------------------
def data_tile():
    path_to_argo = '/net/alpha/exports/sciences/data/ARGO/ARGO/201602-ArgoData'
    path_to_tiles = "/local/tmp/1/herry/tiles/tile%003i.nc"
    path_to_argotiles = '/local/tmp/1/herry/filter/argo_tile%003i.nc'
    return path_to_argo, path_to_tiles, path_to_argotiles


#  ----------------------------------------------------------------------------
def data_argotools():
    path_to_argolog = data_argodb
    path_to_wmo = '/net/alpha/exports/sciences/data/ARGO/ARGO'
    return path_to_wmo, path_to_argolog
