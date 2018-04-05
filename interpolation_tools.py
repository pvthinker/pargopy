# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 14:37:02 2018

@author: herry
"""
import gsw as gsw
import numpy as np
import argotools as argo


def interpolate_profiles(subargodb, wmodic):
    """Interpolate the profiles in subargodb"""

    infos = argo.retrieve_infos_from_tag(subargodb, subargodb['TAG'])

    zref = argo.zref
    n_zref = len(zref)
    n_profiles = len(np.where(subargodb['FLAG'][:] == 0)[0])
    print('nb of valid (flag==0) profiles: %i' % n_profiles)
    CT = np.zeros((n_profiles, n_zref))
    SA = np.zeros((n_profiles, n_zref))
    RHO = np.zeros((n_profiles, n_zref))
    LON = np.zeros((n_profiles,))
    LAT = np.zeros((n_profiles,))
    JULD = np.zeros((n_profiles,))
    TAG = np.zeros((n_profiles,), dtype=int)

    kprof = 0

    wmos = set(infos['WMO'])

    for w in wmos:
        print('interpolate profiles from wmo %i' % w)
        idx = np.where(infos['WMO'] == w)[0]
        print(w, idx)
        iprof = infos['IPROF'][idx]
        dac = argo.dac_from_wmo(wmodic, w)
        data = argo.read_profile(dac, w, header=True, data=True, dataqc=True)
        for l, k in enumerate(iprof):
            if subargodb['FLAG'][idx[l]] == 0:
                temp = data['TEMP'][k, :]
                psal = data['PSAL'][k, :]
                pres = data['PRES'][k, :]
                temp_qc = data['TEMP_QC'][k, :]
                psal_qc = data['PSAL_QC'][k, :]
                pres_qc = data['PRES_QC'][k, :]
                lon = data['LONGITUDE'][k]
                lat = data['LATITUDE'][k]
                Ti, Si, Ri, zCT, zSA, zz, ierr = raw_to_interpolate(temp, psal, pres,
                                         temp_qc, psal_qc, pres_qc,
                                         lon, lat, zref)
                ierr = 0
                if ierr == 0:
                    CT[kprof, :] = Ti
                    SA[kprof, :] = Si
                    RHO[kprof, :] = Ri
                    TAG[kprof] = subargodb['TAG'][idx[l]]
                    LON[kprof] = subargodb['LONGITUDE'][idx[l]]
                    LAT[kprof] = subargodb['LATITUDE'][idx[l]]
                    JULD[kprof] = subargodb['JULD'][idx[l]]
                    kprof += 1
                else:
                    subargodb['FLAG'][idx[l]] = 202
            else:
                pass

    res = {'CT': CT[:kprof, :],
           'SA': SA[:kprof, :],
           'RHO': RHO[:kprof, :],
           'TAG': TAG[:kprof],
           'LONGITUDE': LON[:kprof],
           'LATITUDE': LAT[:kprof],
           'JULD': JULD[:kprof]}

    return res

def raw_to_interpolate(temp, sal, pres, temp_qc, sal_qc, pres_qc, lon, lat, zref):
    """Interpolate in situ data on zref depths
    ierr = 0: no pb
    ierr > 0: pb
    """
    #  print(len(temp_qc))
    CT, SA, z = [], [], []
    klist, ierr = remove_bad_qc(temp, sal, pres, temp_qc, sal_qc, pres_qc)
    if ierr == 0:
        Tis = temp[klist]
        SP = sal[klist]
        p = pres[klist]

        CT, SA, z = insitu_to_absolute(Tis, SP, p, lon, lat, zref)
        Ti, Si = interp_at_zref(CT, SA, z, zref)
        p = gsw.p_from_z(-zref, lat)
        Ri = gsw.rho(Si, Ti, p)
    else:
        Ti, Si, Ri = [], [], []
    return Ti, Si, Ri, CT, SA, z, ierr


def remove_bad_qc(temp, sal, pres, temp_qc, sal_qc, pres_qc):
    """Return the index list of data for which the three qc's are 1
    and the error flag ierr
    ierr = 0 : no pb
    ierr = 1 : too few data in the profile"""
#==============================================================================
#     klist = np.zeros((len(temp_qc))) + n
#     for i, rank in enumerate(temp_qc):
#         if temp_qc[i] == '1':
#             klist[:] = i
#==============================================================================
    klist = [k for k in range(len(temp_qc)) if (temp_qc[k] == '1') and (sal_qc[k] == '1') and (pres_qc[k] == '1')]
    #  print(klist)
    nk = len(klist)
    ierr = 0
    if nk < 5:
        ierr = 1
    #  print(klist)
    return klist, ierr


def insitu_to_absolute(Tis, SP, p, lon, lat, zref):
    """Transform in situ variables to TEOS10 variables"""
    #  SP is in p.s.u.
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, Tis, p)
    z = -gsw.z_from_p(p, lat)
    return(CT, SA, z)


def interp_at_zref(CT, SA, z, zref):
    """Interpolate CT and SA from their native depths z to zref"""
    nbpi, zs, ks = select_depth(zref, z)
    nref = len(zref)
    CTi = np.zeros((nref,), dtype=float)
    SAi = np.zeros((nref,), dtype=float)
    # extrapolation at surface
    k = 0
    nupper = nbpi[0]+nbpi[1]    
    if nupper==1:
        CTi[k] = CT[0]
        SAi[k] = SA[0]
    elif nupper>=2:
        a, b = lincoef(zref[k], (zs[0]+zs[1])[:2])
        CTi[k] = a*CT[0] + b*CT[1]
        SAi[k] = a*SA[0] + b*SA[1]
    else:
        CTi[k] = np.NaN
        SAi[k] = np.NaN
    # level 1
    k = 1
    if (nbpi[k-1] >= 1) and (nbpi[k]+nbpi[k+1] >= 1):
        k0 = ks[k-1][-1]
        k1 = (ks[k]+ks[k+1])[0]
        a, b = lincoef(zref[k], [z[k0], z[k1]])
        CTi[k] = a*CT[k0] + b*CT[k1]
        SAi[k] = a*SA[k0] + b*SA[k1]
    else:
        CTi[k] = np.NaN
        SAi[k] = np.NaN
        
    # interpolation for interior values
    for k in range(2,nref-2):
        if (nbpi[k-2]+nbpi[k-1] >= 2) and (nbpi[k]+nbpi[k+1] >= 2):
            # bilinear interpolation
            k0 = (ks[k-2]+ks[k-1])[-2:]
            k1 = (ks[k]+ks[k+1])[0:2]
            coef, ierr = cubiccoef(zref[k], list(z[k0])+list(z[k1]))
            a, b, c, d = coef
            if ierr==0:
                CTi[k] = a*CT[k0[0]] + b*CT[k0[1]] + c*CT[k1[0]] + d*CT[k1[1]]
                SAi[k] = a*SA[k0[0]] + b*SA[k0[1]] + c*SA[k1[0]] + d*SA[k1[1]]
            else:
                CTi[k] = np.NaN
                SAi[k] = np.NaN
        elif (nbpi[k-2]+nbpi[k-1] >= 2) and (nbpi[k]+nbpi[k+1] == 1):
            k0 = (ks[k-2]+ks[k-1])[-2:]
            k1 = (ks[k]+ks[k+1])[0]
            coef, ierr = paraboliccoef(zref[k], list(z[k0])+[z[k1]])
            a, b, c = coef
            if ierr==0:
                CTi[k] = a*CT[k0[0]] + b*CT[k0[1]] + c*CT[k1]
                SAi[k] = a*SA[k0[0]] + b*SA[k0[1]] + c*SA[k1]
            else:
                CTi[k] = np.NaN
                SAi[k] = np.NaN
            
        elif (nbpi[k-2]+nbpi[k-1] == 1) and (nbpi[k]+nbpi[k+1] >=2):
            k0 = (ks[k-2]+ks[k-1])[-2:]
            k1 = (ks[k]+ks[k+1])[0:2]
            coef, ierr = paraboliccoef(zref[k], [z[k0]]+list(z[k1]))
            a, b, c = coef
            if ierr==0:
                CTi[k] = a*CT[k0] + b*CT[k1[0]] + c*CT[k1[1]]
                SAi[k] = a*SA[k0] + b*SA[k1[0]] + c*SA[k1[1]]
            else:
                CTi[k] = np.NaN
                SAi[k] = np.NaN

            
        elif (nbpi[k-2]+nbpi[k-1] >= 1) and (nbpi[k]+nbpi[k+1] >= 1):
            k0 = (ks[k-2]+ks[k-1])[-1]
            k1 = (ks[k]+ks[k+1])[0]
            a, b = lincoef(zref[k], [z[k0], z[k1]])
            CTi[k] = a*CT[k0] + b*CT[k1]
            SAi[k] = a*SA[k0] + b*SA[k1]            
        else:
            CTi[k] = np.NaN
            SAi[k] = np.NaN

    k = nref-2
    if (nbpi[k-2]+nbpi[k-1] >= 1) and (nbpi[k] >= 1):
        k0 = (ks[k-2]+ks[k-1])[-1]
        k1 = (ks[k])[0]
        a, b = lincoef(zref[k], [z[k0], z[k1]])
        CTi[k] = a*CT[k0] + b*CT[k1]
        SAi[k] = a*SA[k0] + b*SA[k1]
    else:
        CTi[k] = np.NaN
        SAi[k] = np.NaN
    
    # bottom level
    k = nref-1
    nbottom = nbpi[-2]+nbpi[-1]
    if nbottom==1:
        CTi[k] = CT[-1]
        SAi[k] = SA[-1]
    elif nbottom>=2:
        a, b = lincoef(zref[k], (zs[-2]+zs[-1])[-2:])
        CTi[k] = a*CT[-2] + b*CT[-1]
        SAi[k] = a*SA[-2] + b*SA[-1]
    else:
        CTi[k] = np.NaN
        SAi[k] = np.NaN
    return CTi, SAi


def select_depth(zref, z):
    """Return the number of data points we have between successive zref.
    This is used to decide which interpolation is the best: none,
    linear or cubic. For endpoints (zref=0) and bottom(zref=2000) use
    linear extrapolation

    """
    j = 0
    nz = len(z)
    nref = len(zref)
    nbperintervale = np.zeros((nref,), dtype=int)
    zperint = []
    kperint = []
    for k, z0 in enumerate(zref[1:]):
        n = 0
        zs = []
        ks = []
        while (j < nz) and (z[j] < z0):
            n += 1
            zs.append(z[j])
            ks.append(j)
            j += 1
        nbperintervale[k] = n
        zperint.append(zs)
        kperint.append(ks)
    return nbperintervale, zperint, kperint


def lincoef(z0, zs):
    """Weights for linear interpolation at z0 given the two depths in zs"""
    dz = zs[1]-zs[0]
    if abs(dz)<1e-3:
        a, b = 0.5, 0.5
    else:
        a = (zs[1]-z0)/dz
        b = (z0-zs[0])/dz
    return a, b


def paraboliccoef(z0, zs):
    """Weights for cubic interpolation at z0 given the four depths in zs"""
    coef = np.zeros((3,), dtype=float)
    #  print(zs)
    ncoef = len(zs)
    if ncoef < 3:
        print('**** pb with parabolic interp, only %i different zs' % ncoef)
        ierr = 1
    else:
        ierr = 0
    for i in range(ncoef):
        coef[i] = 1.
        for j in range(ncoef):
            if i != j:
                coef[i] *= (z0-zs[j])/(zs[i]-zs[j])
    return coef, ierr


def cubiccoef(z0, zs):
    """Weights for cubic interpolation at z0 given the four depths in zs"""
    coef = np.zeros((4,), dtype=float)
    ncoef = len(set(zs))
    if ncoef < 4:
        print('**** pb with cubic interp, only %i different zs' % ncoef)
        ierr = 1
    else:
        ierr = 0
    for i in range(ncoef):
        coef[i] = 1.
        for j in range(ncoef):
            if i != j:
                coef[i] *= (z0-zs[j])/(zs[i]-zs[j])
    return coef, ierr


def fixqcarray(qc):
    qcm = np.zeros_like(qc, dtype=int)
    shape = np.shape(qcm)
    for j in range(shape[0]):
        for i in range(shape[1]):
            if qc[j,i] == '1':
                qcm[j,i] = 1
    return qcm
