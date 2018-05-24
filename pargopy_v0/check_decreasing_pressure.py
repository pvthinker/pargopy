"""
Created on Mon Mar 12 13:10:24 2018

File used to avoid values error for pressure, temperature and salinity

"""

import param as param
import numpy as np
import argotools as argotools


path_localdata = param.path_to_data
# wmodic = argotools.read_dic('wmodic', path_localdata)
# argodb = argotools.read_dic('argodb', path_localdata)

# maskarraytype = np.ma.core.MaskedArray


def try_to_remove_duplicate_pressure(p):
    idx = [0]+[l+1 for l, x in enumerate(p[1:]) if p[l] < x]
    return idx


# private_list = [1900976]


def check_pressure(p):

    dp = np.diff(p)
    if np.all(dp > 0):
        # print('ok')
        ierr = 0
    else:
        npb = np.sum(np.diff(p) <= 0)
        #  msg = 'dac: %8s / wmo: %8i / iprof: %3i' % (dac, w, l)
        if len(p) < 5:
            print(': only %i p points' % (len(p)))
            ierr = 1
            # print(p)
        else:
            p0 = p.copy()
            idxf = try_to_remove_duplicate_pressure(p)
            p = p[idxf]
            dp = np.diff(p)
            if np.all(dp > 0):
                ierr = 0
                print(': fixed %i pbs' % npb)
                # print(p0, p)
            else:
                if len(p) > 100:
                    # print(p)
                    ierr = 1
                    print(': unfixed')
                else:
                    ierr = 1
                    print(': unfixed [hr profile]')

    return ierr
