# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 13:10:24 2018

@author: herry
"""
# import jdcal
from __future__ import with_statement
import os
import glob
import pickle
import time
import numpy as np
import param as param
import matplotlib.pyplot as plt
import argotools as argotools
tmps1 = time.time()

path_argo = param.path_to_argo
daclist = argotools.daclist
path_localdata = param.path_to_data


def get_all_wmos():
    """Return a dictionnary of all wmo (list of int) with dac (string) as
       keys HAVING a *_prof.nc file

       A few wmo have no *_prof.nc (352 exactly) because ... they have
       actually no profile reported """
    wmodic = {}
    for dac in daclist:
        prfiles = glob.glob('{}/{}/*/*_prof.nc'.format(path_argo, dac))
        wmodic[dac] = [int(f.split('/')[-2]) for f in prfiles]
    return wmodic


#  ----------------------------------------------------------------------------
def get_header_of_all_wmos(wmodic):
    """Get the header of all wmo"""

    n_wmo = argotools.count_wmos(wmodic)
    keys = ['DACID', 'WMO', 'N_PROF', 'N_LEVELS', 'DATE_UPDATE']
    wmostats = {'N_WMO': n_wmo}
    for k in keys:
        wmostats[k] = np.zeros((n_wmo, ), dtype=int)

    iwmo = 0
    for dac in daclist:
        for w in wmodic[dac]:
            output = argotools.read_profile(dac, w, header=True, verbose=False)
            for k in keys:
                wmostats[k][iwmo] = output[k]
            iwmo += 1
            print('{:5,}/{:,} : {} - {}'.format(iwmo, n_wmo, dac, w))

    return wmostats


#  ----------------------------------------------------------------------------
def get_header_of_all_profiles(wmostats):
    """Build argodb from the infos in wmostats

    Once it is created it is more efficient to read it from the disk
    using 'read_argodb()'

    """
    n_profiles = argotools.count_profiles_in_database(wmostats)
    n_wmo = wmostats['N_WMO']
    key_int = ['TAG', 'FLAG']
    key_float = ['LONGITUDE', 'LATITUDE', 'JULD']

    argodb = {}
    for k in key_int:
        argodb[k] = np.zeros((n_profiles,), dtype=int)
    for k in key_float:
        argodb[k] = np.zeros((n_profiles,))

    iprof = 0
    for k in range(n_wmo):
        kdac, wmo = wmostats['DACID'][k], wmostats['WMO'][k]
        dac = daclist[kdac]
        output = argotools.read_profile(dac, wmo, header=True, verbose=False)
        n_prof = output['N_PROF']
        for key in key_float:
            argodb[key][iprof:iprof+n_prof] = output[key]

        for kprof in range(n_prof):
            tag = argotools.get_tag(kdac, wmo, kprof)
            argodb['TAG'][iprof+kprof] = tag

        print(' {:8,}/{:,} : {} - {}'.format(iprof, n_profiles, dac, wmo))
        iprof += n_prof
    return argodb


#  ----------------------------------------------------------------------------
def write_wmodic(wmodic):
    with open('%s/wmodic.pkl' % path_localdata, 'w') as f:
        pickle.dump(wmodic, f)


#  ----------------------------------------------------------------------------
def read_wmodic():
    print('read wmodic.pkl')
    with open('%s/wmodic.pkl' % path_localdata, 'r') as f:
        wmodic = pickle.load(f)
    return wmodic


#  ----------------------------------------------------------------------------
def write_wmstats(wmstats):
    """Read the full wmstats database"""

    with open('%s/wmstats.pkl' % path_localdata, 'w') as f:
        pickle.dump(wmstats, f)


#  ----------------------------------------------------------------------------
def read_wmstats():
    print('read wmostats.pkl')
    with open('%s/wmstats.pkl' % path_localdata, 'r') as f:
        wmstats = pickle.load(f)
    return wmstats


#  ----------------------------------------------------------------------------
def write_argodb(argodb):
    with open('%s/argodb.pkl' % path_localdata, 'w') as f:
        pickle.dump(argodb, f)


#  ----------------------------------------------------------------------------
def read_argodb():
    """Read the full argodb database"""

    print('read argodb.pkl')
    with open('%s/argodb.pkl' % path_localdata, 'r') as f:
        argodb = pickle.load(f)
    return argodb


#  ----------------------------------------------------------------------------
def update_wmodic():
    """Read the full argodb database and update argodb.pkl"""
    wmodic = {}
    new_wmodic = get_all_wmos()
    old_wmodic = read_wmodic()
    for dac in daclist:
        new_wmodic[dac] = set(new_wmodic[dac])
        old_wmodic[dac] = set(old_wmodic[dac])
        wmodic[dac] = new_wmodic[dac].difference(old_wmodic[dac])
    write_wmodic(new_wmodic)


#  ----------------------------------------------------------------------------
if __name__ == '__main__':

    #  Calling the constructor of the class
#==============================================================================
#     if not os.path.isfile('argodb.pkl'):
#         if not os.path.isfile('wmstats.pkl'):
#             if not os.path.isfile('wmodic.pkl'):
#                 wmodic = get_all_wmos()
#                 write_wmodic(wmodic)
#             else:
#                 update_wmodic()
#             wmostats = get_header_of_all_wmos(wmodic)
#             write_wmstats(wmostats)
#         else:
#             pass
#         argodb = get_header_of_all_profiles(wmostats)
#         argofinal = argotools.flag_argodb(argodb, wmostats)
#         write_argodb(argofinal)
#     else:
#         pass
#==============================================================================
    wmodic = get_all_wmos()
    write_wmodic(wmodic)
    wmostats = get_header_of_all_wmos(wmodic)
    write_wmstats(wmostats)
    argodb = get_header_of_all_profiles(wmostats)
    argodb = argotools.flag_argodb(argodb, wmodic)
    write_argodb(argodb)
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
