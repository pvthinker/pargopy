#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 07:32:54 2018

@author: therry

Module contenant la série d'outils utilisée pour l'interpolation des profils
ARGO :
    - Transformation des variables in-situ en variables TEOS-10
    - Interpolation sur les niveaux zref des variables in-situ et TEOS-10

"""

import numpy as np
import pandas as pd
import gsw
from numba import jit

import general_tools as tools
import database as db
import param as param

def interpolate_profiles(subargodb):
    """Interpolate the profiles in subargodb

    :rtype: dic"""

    all_tags_infos = tools.retrieve_infos_from_tag(subargodb.index)
    subargodb['WMO'] = all_tags_infos['WMO']
    subargodb['IDAC'] = all_tags_infos['IDAC']
    subargodb['IPROF'] = all_tags_infos['IPROF']
    
    subargo_to_interp = subargodb[(subargodb['FLAG'] == 0) & (subargodb['STATUS'] == False)]
    nprof_todo = len(subargo_to_interp)
    print('interpolation: #prof: %i, #flag=0/status=False: %i' 
          % (len(subargodb), nprof_todo))
    
    CT = pd.DataFrame(index = subargo_to_interp.index, columns = param.zref)
    SA = pd.DataFrame(index = subargo_to_interp.index, columns = param.zref)
    RHO = pd.DataFrame(index = subargo_to_interp.index, columns = param.zref)
    BVF2 = pd.DataFrame(index = subargo_to_interp.index, columns = param.zref)
    
    set_subargodb = subargo_to_interp.drop_duplicates('WMO')

    counter = 0
    for tag in set_subargodb.index:
        wmo = set_subargodb.loc[tag, 'WMO']
        dac = set_subargodb.loc[tag, 'IDAC']
        part_subargodb = subargo_to_interp[subargo_to_interp['WMO'] == wmo]
        data = db.read_profile(dac, wmo, header=True, data=True, dataqc=True)
        for iprof in part_subargodb['IPROF']:
            print('%i/%i' % (counter, nprof_todo))
            tag_id = part_subargodb.index[part_subargodb['IPROF'] == iprof]
            temp = data['TEMP'][iprof, :]
            psal = data['PSAL'][iprof, :]
            pres = data['PRES'][iprof, :]
            temp_qc = data['TEMP_QC'][iprof, :]
            psal_qc = data['PSAL_QC'][iprof, :]
            pres_qc = data['PRES_QC'][iprof, :]
            lon = data['LONGITUDE'][iprof]
            lat = data['LATITUDE'][iprof]
            Ti, Si, Ri, BVF2i, zCT, zSA, zz, ierr = raw_to_interpolate(temp, psal, pres,
                                                                       temp_qc, psal_qc, pres_qc,
                                                                       lon, lat, param.zref)
            
            del(temp, psal, pres, temp_qc, psal_qc, pres_qc, lon, lat)
            
            ierr = 0
            if len(Ti) == 0:
                ierr = 1
            # Checking if Ti, Si, Ri are full of NaN
            checker = []
            checker2 = []
            for i, T in enumerate(Ti):
                checker.append(np.isnan(T))
                checker2.append(np.isnan(Si[i]))
            if (False in checker) and (False in checker2):
                ierr = 0
            else:
                ierr = 2
            if ierr == 0:
                CT.loc[tag_id, :] = Ti
                SA.loc[tag_id, :] = Si
                RHO.loc[tag_id, :] = Ri
                BVF2.loc[tag_id, :] = BVF2i
            else:
                subargodb.loc[tag_id, 'FLAG'] = 202
            counter += 1
        else:
            pass
        
    # CT = CT.dropna(how='all')
    # SA = SA.dropna(how='all')
    # RHO = RHO.dropna(how='all')
    # BVF2 = BVF2.dropna(how='all')
    
    res = {'CT': CT,
           'SA': SA,
           'RHO': RHO,
           'BVF2': BVF2}
    subargodb['STATUS'] = True
    print('Interpolation ended')
    return res, subargodb


def raw_to_interpolate(temp, sal, pres, temp_qc, sal_qc, pres_qc, lon, lat, zref):
    """Interpolate in situ data on zref depths

    ierr = 0: no pb

    ierr > 0: pb

    :rtype: float, float, float, float, float, float, float, int
    """

    klist, ierr = remove_bad_qc(temp, sal, pres, temp_qc, sal_qc, pres_qc)
    if ierr == 0:
        Tis = temp[klist]
        SP = sal[klist]
        p = pres[klist]

        CT, SA, z = insitu_to_absolute(Tis, SP, p, lon, lat, zref)
        Ti, Si, dTidz, dSidz = interp_at_zref(CT, SA, z, zref)
        pi = gsw.p_from_z(-zref, lat)
        #Ri = gsw.rho(Si, Ti, pi)
        Ri, alpha, beta = gsw.rho_alpha_beta(Si, Ti, pi)
        g = gsw.grav(lat, pi)
        BVF2i = g*(beta*dSidz-alpha*dTidz)
        # Nf = np.sqrt(N2)/gsw.f(lat)
    else:
        Ti, Si, Ri, BVF2i, CT, SA, z= [], [], [], [], [], [], []
    return Ti, Si, Ri, BVF2i, CT, SA, z, ierr


def remove_bad_qc(temp, sal, pres, temp_qc, sal_qc, pres_qc):
    """Return the index list of data for which the three qc's are 1
    and the error flag ierr

    ierr = 0 : no pb

    ierr = 1 : too few data in the profile

    :rtype: list, int"""
    maskarraytype = np.ma.core.MaskedArray
    keys = [temp, sal, pres, temp_qc, sal_qc, pres_qc]
    for key in keys:
        if key.any == maskarraytype:
            key = key.compressed()
    klist = [k for k in range(len(pres)) if (temp_qc[k] == '1') and (
        sal_qc[k] == '1') and (pres_qc[k] == '1')]
    ierr = 0
    p = pres[klist]
    ierr = check_pressure(p)

    return klist, ierr


def insitu_to_absolute(Tis, SP, p, lon, lat, zref):
    """Transform in situ variables to TEOS10 variables

    :rtype: float, float, float"""
    #  SP is in p.s.u.
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, Tis, p)
    z = -gsw.z_from_p(p, lat)
    return(CT, SA, z)


def interp_at_zref(CT, SA, z, zref):
    """Interpolate CT, SA, dCT/dz and dSA/dz from their native depths z to
    zref

    Method: we use piecewise Lagrange polynomial interpolation

    For each zref[k], we select a list of z[j] that are close to
    zref[k], imposing to have z[j] that are above and below zref[k]
    (except near the boundaries)

    If only two z[j] are found then the result is a linear interpolation

    If n z[j] are found then the result is a n-th order interpolation.

    For interior points we may go up to 6-th order

    For the surface level (zref==0), we do extrapolation

    For the bottom level (zref=2000), we do either extrapolation or
    interpolation if data deeper than 2000 are available.

    :rtype: float, float, float, float

    """

    nref = len(zref)
    CTi = np.zeros((nref,), dtype=float)
    SAi = np.zeros((nref,), dtype=float)
    dCTdzi = np.zeros((nref,), dtype=float)
    dSAdzi = np.zeros((nref,), dtype=float)

    nbpi, ks = select_depth(zref, z)
    nupper = np.zeros((nref,), dtype=int)
    nlower = np.zeros((nref,), dtype=int)

    # count the number of data that are lower and upper than zref[k]
    for k in range(nref):
        if k > 0:
            nlower[k] += nbpi[k-1]
        if k > 1:
            nlower[k] += nbpi[k-2]
        if k < nref:
            nupper[k] += nbpi[k]
        if k < nref-1:
            nupper[k] += nbpi[k+1]

    # for each zref, form the list of z[j] used for the interpolation
    # if the list has at least two elements (a linear interpolation is possible)
    # then do it, otherwise, skip that depth
    for k in range(nref):
        idx = []
        if k == 0:
            if nupper[k] >= 2:
                idx = (ks[0]+ks[1]+ks[2])[:3]
        elif k == 1:
            if (nlower[k] >= 1) and (nupper[k] >= 1):
                idx = ks[0][:-3]+(ks[1]+ks[2])[:3]
        elif k == (nref-1):
            if (nlower[k]+nupper[k]) >= 2:
                idx = (ks[k-2]+ks[k-1])[-3:]+ks[k][:3]
        elif k == (nref-2):
            if (nlower[k] >= 1) and (nupper[k] >= 1):
                idx = (ks[k-2]+ks[k-1])[-3:]+(ks[k]+ks[k+1])[:3]
        else:
            if (nlower[k] >= 1) and (nupper[k] >= 1):
                idx = (ks[k-2]+ks[k-1])[-3:]+(ks[k]+ks[k+1])[:3]

        if len(idx) >= 2:
            cs, ds = lagrangepoly(zref[k], z[idx])
            # the meaning of the weights computed by lagrangepoly should
            # be clear in the code below
            #
            # cs[i] (resp. ds[i]) is the weight to apply on CT[idx[i]]
            # sitting at z[idx[i]] to compute CT (resp. dCT/dz) at zref[k]
            #
            CTi[k] = np.sum(cs*CT[idx])
            SAi[k] = np.sum(cs*SA[idx])
            dCTdzi[k] = np.sum(ds*CT[idx])
            dSAdzi[k] = np.sum(ds*SA[idx])
        else:
            CTi[k] = np.nan
            SAi[k] = np.nan
            dCTdzi[k] = np.nan
            dSAdzi[k] = np.nan

    return CTi, SAi, dCTdzi, dSAdzi


@jit
def select_depth(zref, z):
    """Return the number of data points we have between successive zref.

    for each intervale k, we select the z_j such that

    zref[k] <= z_j < zref[k+1], for k=0 .. nref-2

    zref[nref-1] <= z_j < zextra, for k=nref-1

    and return

    nbperintervale[k] = number of z_j

    kperint[k] = list of j's


    with zextra = 2*zref[-1] - zref[-2]

    :rtype: int, list

    """
    nz = len(z)
    nref = len(zref)
    zextra = 2*zref[-1]-zref[-2]
    zrefextended = list(zref)+[zextra]
    nbperintervale = np.zeros((nref,), dtype=int)
    kperint = []
    zprev = -1.
    j = 0
    for k, z0 in enumerate(zrefextended[1:]):
        n = 0
        ks = []
        while (j < nz) and (z[j] < z0):
            # for a few profiles it may happens that two consecutive
            # data sit at the same depth this causes a division by
            # zero in the interpolation routine.  Here we fix this by
            # simply skipping depths that are already used.
            if z[j] > zprev:
                n += 1
                ks.append(j)
            zprev = z[j]
            j += 1
        nbperintervale[k] = n
        kperint.append(ks)
    return nbperintervale, kperint


@jit
def lagrangepoly(x0, xi):
    """Weights for polynomial interpolation at x0 given a list of xi
    return both the weights for function (cs) and its first derivative
    (ds)

    Example:
    lagrangepoly(0.25, [0, 1])
    >>> [0.75, 0.25,], [1, -1]

    :rtype: float, float

    """
    xi = np.asarray(xi)
    ncoef = len(xi)
    cs = np.ones((ncoef,))
    ds = np.zeros((ncoef,))

    denom = np.zeros((ncoef, ncoef))
    for i in range(ncoef):
        for j in range(ncoef):
            if i != j:
                dx = xi[i]-xi[j]
                if dx == 0:
                    # should not happen because select_depth removes
                    # duplicate depths
                    #  raise ValueError('division by zero in lagrangepoly')
                    print('WARNING, division by zero in lagrangepoly')
                else:
                    denom[i, j] = 1./dx

    for i in range(ncoef):
        for j in range(ncoef):
            if i != j:
                cff = 1.
                cs[i] *= (x0-xi[j])*denom[i, j]
                for k in range(ncoef):
                    if (k != i) and (k != j):
                        cff *= (x0-xi[k])*denom[i, k]
                ds[i] += cff*denom[i, j]
    return cs, ds


@jit
def try_to_remove_duplicate_pressure(p):
    """
    :param p: list of pressures
    
    Fonction utilisée pour trouver d'éventuelles pressions dupliquées
    
    :rtype: list
    """
    idx = [0]+[l+1 for l, x in enumerate(p[1:]) if p[l] < x]
    return idx


@jit
def check_pressure(p):
    """
    :param p: list of pressures
    
    Fonction utilisée pour vérifier les valeurs de pression et les flagger comme
    mauvaises si elles ne conviennent pas
    
    :rtype: int
    """

    
    dp = np.diff(p)
    if np.all(dp > 0):
        ierr = 0
    else:
        npb = np.sum(np.diff(p) <= 0)
        if len(p) < 5:
            print(': only %i p points' % (len(p)))
            ierr = 1
            # print(p)
        else:
            idxf = try_to_remove_duplicate_pressure(p)
            p = p[idxf]
            dp = np.diff(p)
            if np.all(dp > 0):
                ierr = 0
                print(': fixed %i pbs' % npb)
                # print(p0, p)
            else:
                if len(p) > 100:
                    ierr = 1
                    print(': unfixed')
                else:
                    ierr = 1
                    print(': unfixed [hr profile]')

    return ierr
