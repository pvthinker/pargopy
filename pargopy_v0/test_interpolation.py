import gsw as gsw
import numpy as np
import argotools as argotools
import matplotlib.pyplot as plt

from interpolation_tools import *

plt.ion()

itile = 50
wmodic = argotools.read_wmodic()
subargodb = argotools.read_argo_filter(itile)

infos = argotools.retrieve_infos_from_tag(subargodb['TAG'])
wmos = list(set(infos['WMO']))

zref = argotools.zref
n_zref = len(zref)

nwmos = len(wmos)
iwmo = np.random.permutation(np.arange(nwmos))

plt.clf()
for kwmo in iwmo[:5]:

    w = wmos[kwmo]
    dac = argotools.dac_from_wmo(wmodic, w)
    data = argotools.read_profile(dac, w, header=True, data=True, dataqc=True)
    nprof = len(data['JULD'])

    iprof = np.random.permutation(np.arange(nprof))
    for k in iprof[:5]:

        temp = data['TEMP'][k, :]
        psal = data['PSAL'][k, :]
        pres = data['PRES'][k, :]
        temp_qc = data['TEMP_QC'][k, :]
        psal_qc = data['PSAL_QC'][k, :]
        pres_qc = data['PRES_QC'][k, :]
        lon = data['LONGITUDE'][k]
        lat = data['LATITUDE'][k]

        Ti, Si, Ri, BVF2, zCT, zSA, zz, ierr = raw_to_interpolate(temp, psal, pres,
                                                                  temp_qc, psal_qc, pres_qc,
                                                                  lon, lat, zref)

        if ierr == 0:
            plt.subplot(121)
            #plt.plot(temp, -pres, '.')
            plt.plot(Ti, -zref)
            # plt.plot(zCT.data, -pres.data)
            ax = plt.gca()
            ax.set_ylim([-2050, 0])

            plt.subplot(122)
            plt.semilogx(BVF2, -zref)
            plt.axis([1e-8, 1e-2, -2050, 0])
