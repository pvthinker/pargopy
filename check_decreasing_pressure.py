import pickle
import research_tools as research
import argotools as argotools
import param as param
import numpy as np

wmodic = argotools.read_wmodic()
argodb = argotools.read_argodb()

nprof = len(argodb['JULD'])

infos = argotools.retrieve_infos_from_tag(argodb, argodb['TAG'])

wmos = list(set(infos['WMO']))

maskarraytype = np.ma.core.MaskedArray

def try_to_remove_duplicate_pressure(p):
    idx = [0]+[l+1 for l, x in enumerate(p[1:]) if p[l]<x]
    return idx

private_list = [1900976]

for w in wmos:
    idx = np.where(infos['WMO'] == w)[0]
    iprof = infos['IPROF'][idx]
    dac = argotools.dac_from_wmo(wmodic, w)
    data = argotools.read_profile(dac, w, header=True,
                                  data=True, dataqc=True, verbose=False)
    nlevs = data['N_LEVELS']
    for l, k in enumerate(iprof):
        if (argodb['FLAG'][idx[l]] == 0):
            p = data['PRES'][l, :]
            if type(p) == maskarraytype:
                p = p.compressed()
            else:
                pass
                # print('dac: %s / wmo: %i / iprof: %i => not masked' % (dac, w, l))
                
            dp = np.diff(p)
            if np.all(dp > 0):
                # print('ok')
                pass
            else:
                npb = np.sum(np.diff(p)<=0)
                msg = 'dac: %8s / wmo: %8i / iprof: %3i' % (dac, w, l)
                if len(p) < 5:
                    print(msg+': only %i p points' % (len(p)))
                    # print(p)
                else:
                    p0 = p.copy()
                    idxf = try_to_remove_duplicate_pressure(p)
                    p = p[idxf]
                    dp = np.diff(p)
                    if np.all(dp > 0):
                        print(msg + ': fixed %i pbs' % npb)
                        # print(p0, p)
                    else:
                        if nlevs > 100:
                            #print(p)
                            print(msg + ': unfixed')
                        else:
                            print(msg + ': unfixed [hr profile]')
                        # print(p)

