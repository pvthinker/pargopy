# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 13:10:24 2018

File creating the summary of ARGO used to generate the atlas

"""
# import jdcal
from __future__ import with_statement
import glob
import time
import numpy as np
import param as param
#  import matplotlib.pyplot as plt
import argotools as argotools
tmps1 = time.time()

path_argo = param.path_to_argo
daclist = ['aoml', 'bodc', 'coriolis', 'csio',
           'csiro', 'incois', 'jma', 'kma',
           'kordi', 'meds', 'nmdis']
path_localdata = param.path_to_data


def get_all_wmos():
    """
    Return a dictionnary of all wmo (list of int) with dac (string) as
    keys HAVING a \*_prof.nc file

    A few wmo have no \*_prof.nc (352 exactly) because ... they have
    actually no profile reported 
    
    :rtype: dic
    """
    wmodic = {}
    for dac in daclist:
        prfiles = glob.glob('{}/{}/*/*_prof.nc'.format(path_argo, dac))
        wmodic[dac] = [int(f.split('/')[-2]) for f in prfiles]
    return wmodic


#  ----------------------------------------------------------------------------
def get_header_of_all_wmos(wmodic):
    """
    Get the header of all wmo
    
    :rtype: dic
    """

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
    
    :rtype: dic

    """
    n_profiles = argotools.count_profiles_in_database(wmostats)
    n_wmo = wmostats['N_WMO']
    key_int = ['TAG', 'FLAG']
    key_float = ['LONGITUDE', 'LATITUDE', 'JULD']
    key_char = ['DATA_MODE']

    argodb = {}
    for k in key_int:
        argodb[k] = np.zeros((n_profiles,), dtype=int)
    for k in key_float:
        argodb[k] = np.zeros((n_profiles,))
    for k in key_char:
        argodb[k] = np.zeros((n_profiles,), dtype='c')

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

        for i, key in enumerate(key_char):
            argodb[key][iprof:iprof+n_prof] = output[key]

        print(' {:8,}/{:,} : {} - {}'.format(iprof, n_profiles, dac, wmo))
        iprof += n_prof
    return argodb


#  ----------------------------------------------------------------------------
def update_wmodic():
    """
    Read the full argodb database and update argodb.pkl
    
    :rtype: None
    """
    wmodic = {}
    new_wmodic = get_all_wmos()
    #  old_wmodic = argotools.read_wmodic()
    old_wmodic = argotools.read_dic('wmodic', path_localdata)
    for dac in daclist:
        new_wmodic[dac] = set(new_wmodic[dac])
        old_wmodic[dac] = set(old_wmodic[dac])
        wmodic[dac] = new_wmodic[dac].difference(old_wmodic[dac])
    argotools.write_dic('wmodic', new_wmodic, path_localdata)


#  ----------------------------------------------------------------------------
def propagate_flag_backward(argodb, subargodb, verbose=True):
    """
    Update argodb FLAG using subargodb
    :rtype: None
    """

    for k, tag in enumerate(subargodb['TAG']):
        idx = np.where(argodb['TAG'] == tag)[0][0]
        prev = argodb['FLAG'][idx]
        new = subargodb['FLAG'][k]
        if prev == new:
            pass
        else:
            if verbose:
                print('tag %i flag changed from %i to %i' % (tag, prev, new))
            argodb['FLAG'][idx] = subargodb['FLAG'][k]
    print('Going to rewrite argodb')
    argotools.write_dic('argodb', argodb, path_localdata)
    print('Argodb rewrited')


#  ----------------------------------------------------------------------------
def main():
    """Main function of argodb.py"""
    wmodic = get_all_wmos()
    #  write_wmodic(wmodic)
    argotools.write_dic('wmodic', wmodic, path_localdata)
    wmostats = get_header_of_all_wmos(wmodic)
    #  write_wmstats(wmostats)
    argotools.write_dic('wmostats', wmostats, path_localdata)
    argodb = get_header_of_all_profiles(wmostats)
    argodb = argotools.flag_argodb(argodb, wmodic)
    #  write_argodb(argodb)
    argotools.write_dic('argodb', argodb, path_localdata)


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
    main()
    tmps2 = time.time() - tmps1
    print("Temps d'execution = %f" % tmps2)
