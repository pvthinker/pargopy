# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 11:01:25 2018

@author: herry

File defining the pathes to use this program

"""
#  location is used to know where you are working
#  if you're working on datarmor, location takes the value 'DATARMOR'
#  if you're working on herry-s local, location takes the value 'HERRY'

#  location = 'HERRY'
location = 'DATARMOR_TMP'
#  location = 'DATARMOR_FINAL'

if location == 'DATARMOR_TMP':
    # path_to_argo = '/home/oo26/coriolis/co05/co0508/dac'
    path_to_argo = '/datawork/fsi2/coriolis-s/public/co05/co0508/gdac/dac'

    path_to_pargopy = '/home2/datahome/therry/pargopy/pargopy_v0/'

    path_to_data = '/home2/datawork/therry/data'

    path_to_filter = '/home2/datawork/therry/filter'

    path_to_stats = '/home2/datawork/therry/tmp/stats'

    path_to_tiles = '/home2/datawork/therry/tmp/tiles'

    path_to_atlas = '/home2/datawork/therry/tmp/atlas'

if location == 'DATARMOR_FINAL':
    # path_to_argo = '/home/oo26/coriolis/co05/co0508/dac'
    path_to_argo = '/datawork/fsi2/coriolis-s/public/co05/co0508/gdac/dac'

    path_to_pargopy = '/home2/datahome/therry/pargopy/pargopy_v0/'

    path_to_data = '/home2/datawork/therry/data'

    path_to_filter = '/home2/datawork/therry/filter'

    path_to_stats = '/home2/datawork/therry/final/stats'

    path_to_tiles = '/home2/datawork/therry/final/tiles'

    path_to_atlas = '/home2/datawork/therry/final/atlas'

elif location == 'HERRY':
    path_to_argo = '/net/alpha/exports/sciences/data/ARGO/ARGO/201602-ArgoData'

    path_to_pargopy = '/home/herry/Documents/pargopy/pargopy_v0/'

    path_to_data = '/local/tmp/1/herry/pargopy/data'

    path_to_filter = '/local/tmp/1/herry/pargopy/filter'

    path_to_stats = '/local/tmp/1/herry/pargopy/stats'

    path_to_tiles = '/local/tmp/1/herry/pargopy/tiles'

    path_to_atlas = '/local/tmp/1/herry/pargopy/atlas'


elif location == 'roullet_lops':
    path_to_argo = '/net/alpha/exports/sciences/data/ARGO/ARGO/201602-ArgoData'

    path_to_pargopy = '/home/roullet/dev/pargopy'

    workdir = '/net/libra/local/tmp/1/roullet/data/ARGO/pargopy'

    path_to_data = '%s/data' % workdir

    path_to_filter = '%s/filter' % workdir

    path_to_stats = '%s/stats' % workdir

    path_to_tiles = '%s/tiles' % workdir

    path_to_atlas = '%s/atlas' % workdir
