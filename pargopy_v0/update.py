# import pickle
import argotools as at
import glob
import pandas as pd
import os

# f = open('argodb')
# a = pickle.load(f)
# f.close()

# tag = a['TAG']

# infos = at.retrieve_infos_from_tag(a, tag)
# wmos = infos['WMO']

path_argo = "/home/ref-argo/gdac/dac"

daclist = ['aoml', 'bodc', 'coriolis', 'csio',
           'csiro', 'incois', 'jma', 'kma',
           'kordi', 'meds', 'nmdis']


def get_wmo_infos(path_gdac):
    w = []
    d = []
    n = []
    u = []
    for dac in daclist:  # ['coriolis']
        prfiles = glob.glob('{}/{}/*/*_prof.nc'.format(path_gdac, dac))
        wmos = [int(f.split('/')[-2]) for f in prfiles]
        for wmo in wmos:
            output = at.read_profile(dac, wmo, shortheader=True, verbose=True,
                                     path=path_gdac)
            w += [wmo]
            d += [dac]
            n += [output['N_PROF']]
            u += [output['DATE_UPDATE']]

    return pd.DataFrame({'DAC': d, 'N_PROF': n, 'DATA_UPDATE': u}, index=w)


if os.path.isfile('oldwmo'):
    old = pd.read_pickle('oldwmo')
else:
    old = get_wmo_infos(at.param.path_to_argo)
    old.to_pickle('oldwmo')

if os.path.isfile('newwmo'):
    new = pd.read_pickle('newwmo')
else:
    new = get_wmo_infos(path_argo)
    new.to_pickle('newwmo')

verbose = False
# wmos counters
nchanged = 0
nidentical = 0
nnew = 0
# profiles counters
pchanged = 0
pidentical = 0
pnew = 0
print('Compare the two gdac repositories')
for w in new.index:
    if w in old.index:
        changed = any(new.loc[w, :] != old.loc[w, :])
        if changed:
            nchanged += 1
            if new.loc[w, 'N_PROF'] == old.loc[w, 'N_PROF']:
                if verbose:
                    print('%i has been updated' % w)
                pchanged += new.loc[w, 'N_PROF']
            else:
                if verbose:
                    print('%i has more profiles : %i => %i' %
                          (w, old.loc[w, 'N_PROF'], new.loc[w, 'N_PROF']))
                pnew += new.loc[w, 'N_PROF']-old.loc[w, 'N_PROF']
                pidentical += old.loc[w, 'N_PROF']
        else:
            nidentical += 1
            pidentical += new.loc[w, 'N_PROF']
    else:
        if verbose:
            print('%i is new' % w)
        nnew += 1
        pnew += new.loc[w, 'N_PROF']

print('-'*60)
print('Identical wmos : %i' % nidentical)
print('Changed   wmos : %i' % nchanged)
print('New       wmos : %i' % nnew)
print('-'*60)
print('Identical profiles : %i' % pidentical)
print('Changed   profiles : %i' % pchanged)
print('New       profiles : %i' % pnew)
print('-'*60)

diff_wmos = pd.merge(new, old, how='inner')

# new_wmos = []
# for dac in daclist:
#     new_wmos += wmodic[dac]

# w0 = set(wmos)
# w1 = set(new_wmos)

# wold = w1.intersection(w0)
# wnew = w1.difference(w0)

# # check in old wmos who has new profiles
