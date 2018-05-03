# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 13:10:24 2018

@author: herry

Tools used by all the python files to generate the atlas

"""

from __future__ import with_statement
import jdcal
import os
import glob
import pickle
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import param as param
import melted_functions as melted

path_argo = param.path_to_argo
path_filter = param.path_to_filter
path_localdata = param.path_to_data

daclist = ['aoml', 'bodc', 'coriolis', 'csio',
           'csiro', 'incois', 'jma', 'kma',
           'kordi', 'meds', 'nmdis']

key_header = ['LATITUDE', 'LONGITUDE', 'JULD']
key_headerqc = ['POSITION_QC', 'JULD_QC']

# home-made QC telling whether a profile has 'TEMP', 'PSAL' and 'PRES'
# present. A few profiles don't have 'PSAL'
key_variableqc = ['TSP_QC']

key_data = ['TEMP', 'PSAL', 'PRES']
key_dataqc = ['TEMP_QC', 'PSAL_QC', 'PRES_QC']

zref = np.array([0., 10., 20., 30., 40., 50., 60., 70., 80., 90.,
                 100., 110., 120., 130., 140., 150., 160., 170.,
                 180., 190., 200., 220., 240., 260., 280, 300.,
                 320., 340., 360., 380., 400., 450., 500., 550.,
                 600., 650., 700., 750., 800., 850., 900., 950.,
                 1000., 1050., 1100., 1150., 1200., 1250., 1300.,
                 1350., 1400., 1450., 1500., 1550., 1600., 1650.,
                 1700., 1750., 1800., 1850., 1900., 1950.,
                 2000.])


def dac_from_wmo(wmodic, wmo):
    """
    Retrieve the dac of a wmo
    :rtype: list of dac
    """

    dac = ''
    for d in daclist:
        if wmo in wmodic[d]:
            dac = d
    return dac


def count_wmos(wmodic):
    """
    Count the total number of wmo in the Argo database base
    :rtype: list of wmo
    """

    nwmos = 0
    for dac in daclist:
        nwmos += len(wmodic[dac])
    print('number of wmo in Argo database: %i' % nwmos)
    return nwmos


def count_profiles_in_database(wmostats):
    """
    Count the total number of profiles in database
    :rtype: int
    """

    nbprofiles = 0
    for nbpr in wmostats['N_PROF']:
        nbprofiles += nbpr
    print('number of profiles in Argo database: %i' % nbprofiles)
    return nbprofiles


def get_profile_file_path(dac, wmo):
    """
    Return the file path to the \*_prof.nc data file
    
    :rtype: string
    """

    if type(dac) is int:
        dac = daclist[dac]
    filename = '%s/%s/%i/%i_prof.nc' % (path_argo, dac, wmo, wmo)
    return filename


def read_profile(dac, wmo, iprof=None,
                 header=False, data=False,
                 headerqc=False, dataqc=False,
                 verbose=True):
    """Basic driver to read the \*_prof.nc data file

    The output is a dictionnary of vectors
    - read one or all profiles read the header (lat, lon, juld) or not
    - read the data or not always return IDAC, WMO, N_PROF, N_LEVELS
    - and DATA_UPDATE (all 5 are int)
    
    :rtype: dic
    """
    if type(dac) is int:
        dac = daclist[dac]

    filename = get_profile_file_path(dac, wmo)

    if verbose:
        print('/'.join(filename.split('/')[-3:]))

    output = {}

    required_keys = set(['TEMP', 'PSAL', 'PRES'])

    if (os.path.isfile(filename)):
        with Dataset(filename, "r", format="NETCDF4") as f:
            output['DACID'] = daclist.index(dac)
            output['WMO'] = wmo
            output['N_PROF'] = len(f.dimensions['N_PROF'])
            output['N_LEVELS'] = len(f.dimensions['N_LEVELS'])
            # DATE_UPDATE is an array of 14 characters in the *_prof.nc
            # we transform it into an int
            # YYYYMMDDhhmmss
            output['DATE_UPDATE'] = ''.join(f.variables['DATE_UPDATE'][:])

            keyvar = set(f.variables.keys())

            if required_keys.issubset(keyvar):
                output['TSP_QC'] = '1'
            else:
                output['TSP_QC'] = '2'

            if header or headerqc or data or dataqc:
                if iprof is None:
                    idx = range(output['N_PROF'])
                    output['IPROF'] = np.arange(output['N_PROF'])
                else:
                    idx = iprof
                    output['IPROF'] = iprof

            if header:
                for key in key_header:
                    output[key] = f.variables[key][idx]
                    output['DATA_MODE'] = np.asarray([c for c in f.variables['DATA_MODE'][idx]])

            if headerqc:
                for key in key_headerqc:
                    output[key] = f.variables[key][idx]

            if data:
                for key in key_data:
                    if output['TSP_QC'] == '1':
                        output[key] = f.variables[key][idx, :]
                    else:
                        output[key] = np.NaN+np.zeros(
                                (output['N_PROF'], output['N_LEVELS']))

            if dataqc:
                for key in key_dataqc:
                    if output['TSP_QC'] == '1':
                        output[key] = f.variables[key][idx]
                    else:
                        output[key] = np.zeros(
                                (output['N_PROF'], output['N_LEVELS']), dtype=str)

    return output


def flag_argodb(argodb, wmodic):
    """
    Add the flag to argodb
    
    :rtype: dic
    """

    infos = retrieve_infos_from_tag(argodb, argodb['TAG'])
    wmos = set(infos['WMO'])
    n_profiles = len(argodb['JULD'])
    argodb['FLAG'] = np.zeros((n_profiles,), dtype=int)
    idx = []
    for w in wmos:
        print('add flag to wmo %i' % w)
        idx = np.where(infos['WMO'] == w)[0]
        iprof = infos['IPROF'][idx]
        dac = dac_from_wmo(wmodic, w)
        res = read_profile(dac, w, headerqc=True)
        # print('-'*80)
        # print(res['TEMP_QC'].data)
        # print(res['POSITION_QC'].data)
        for k in iprof:
            flag = 0
            if res['POSITION_QC'][k] != '1':
                flag = 101
            if res['JULD_QC'][k] != '1':
                flag = 102
            if res['TSP_QC'] != '1':
                flag = 103
            argodb['FLAG'][idx[k]] = flag
    return argodb


def fix_flag_latlonf(argodb):
    """Set a flag error for profiles having bad masked positions

    Bad masked position yield value of 99999. This fix only concerns a
    few profiles from dac='jma'
    
    :rtype: None
    """
    idx = np.where(argodb['LONGITUDE'] == 99999.)[0]
    if len(idx) > 0:
        argodb['FLAG'][idx] = 104
    # no return needed because argodb is mutable


def get_tag(kdac, wmo, kprof):
    """Compute the tag number of a profile

    The inverse of get_tag() is retrieve_infos_from_tag()
    
    :rtype: int
    """
    if kprof > 1000:
        raise ValueError("kprof > 1000, the tag may be wrong")

    return (kdac*10000000+wmo)*1000+kprof


def retrieve_infos_from_tag(argodb, tag):
    """Retrieve idac, wmo and iprof from tag (array of int)

    It is the inverse of get_tag()
    
    :rtype: dic
    """

    iprof = tag % 1000
    tag = (tag-iprof) // 1000
    wmo = tag % 10000000
    tag = (tag-wmo) // 10000000
    idac = tag
    output = {'IDAC': idac, 'WMO': wmo, 'IPROF': iprof}
    return output


def get_datamode(data):
    """Return the data mode of the profile
    
    :rtype: np.array
    """
    mode = np.zeros(len(data), dtype=int)
    for i, d in enumerate(data):
        if d == 'R':
            mode[i] = 0
        elif d == 'A':
            mode[i] = 1
        elif d == 'D':
            mode[i] = 2
        else:
            raise ValueError("No data mode given")
    return mode


def plot_location_profiles(argodb):
    """Plot a scatter plot of profiles in argodb

    argodb can be the full database or any subset
    
    :rtype: None
    """
    idx = np.where(argodb['FLAG'][:] == 0)[0]
    lon = argodb['LONGITUDE'][idx]
    lat = argodb['LATITUDE'][idx]
    plt.figure()
    plt.plot(lon, lat, '.')
    plt.show()


def plot_wmo_data(dac, wmo):
    """Plot raw 'TEMP' data for dac, wmo
    
    :rtype: None"""

    data = read_profile(dac, wmo, data=True)

    varname = 'TEMP'
    plt.figure()
    plt.imshow(data[varname], interpolation='nearest', origin='lower')
    plt.axis('tight')
    plt.title('%s / wmo=%i / %s' % (dac, wmo, varname))
    plt.colorbar()
    plt.show()


def plot_wmos_stats(wmostats):
    """Plot the histogram of number of profiles per number of levels
    
    :rtype: None"""

    plt.figure()
    plt.hist(wmostats['N_LEVELS'],
             weights=wmostats['N_PROF'],
             bins=np.arange(0, 2000, 10))
    plt.xlabel('nb of levels')
    plt.ylabel('nb of profiles')
    plt.show()


#  ----------------------------------------------------------------------------
def read_wmstats():
    """Read the full wmstats database
    
    :rtype: dic"""
    print('read wmostats.pkl')
    with open('%s/wmstats.pkl' % path_localdata, 'r') as f:
        wmstats = pickle.load(f)
    return wmstats


#  ----------------------------------------------------------------------------
def read_argodb():
    """Read the full argodb database
    
    :rtype: dic"""

    print('read argodb.pkl')
    with open('%s/argodb.pkl' % path_localdata, 'rb') as f:
        argodb = pickle.load(f)
    return argodb


#  ----------------------------------------------------------------------------
def read_wmodic():
    """Read the database containing the list of all the wmos
    
    :rtype: dic"""
    print('read wmodic.pkl')
    with open('%s/wmodic.pkl' % path_localdata, 'r') as f:
        wmodic = pickle.load(f)
    return wmodic


#  ----------------------------------------------------------------------------
def read_argo_filter(i):
    """Read the files containing parts of argodb chosen with lat/lon filters
    
    :rtype: dic"""
    print('read argodic%003i.pkl' % i)
    with open('%s/argodic%003i.pkl' % (path_filter, i), 'rb') as f:
        argodic = pickle.load(f)
    return argodic


#  ----------------------------------------------------------------------------
def tile_definition():
    """Creating the variables
    
    :rtype: float, float, int, int, float, float"""

    minlon = -180.
    maxlon = 180.
    latmin = -80.
    maxlat = 80.

    nlon = 20
    nlat = 15

    minmargin = 1.

    deltalon = maxlon-minlon
    deltalat = maxlat-latmin

    lon = minlon + np.arange(nlon+1) * deltalon/nlon
    lat = latmin + np.arange(nlat+1) * deltalat/nlat

    # distribute the latitudes so that their differences
    # vary in cos(lat)
    # do it via an iterative method
    for k in range(5):
        latm = 0.5*(lat[1:]+lat[:-1])
        dlat = np.diff(lat) * np.cos(latm*np.pi/180)
        dlat = dlat*0 + np.mean(dlat)
        dlat = dlat / np.cos(latm*np.pi/180)
        dlat = dlat/sum(dlat)*deltalat
        lat[1:] = latmin + np.cumsum(dlat)

    marginlat = minmargin / np.cos(latm*np.pi/180)
    marginlon = 2

    return lat, lon, nlat, nlon, marginlat, marginlon


#  ----------------------------------------------------------------------------
def conversion_juld_gregd(juld):
    """Method converting julian day into gregorian day
    
    :rtype: list of int"""
    #  lats, lons, juld = self.reading_variables()
    #  2433282.5000000 corresponds to the Argo origin date
    gregday = jdcal.jd2jcal(2433282.500000, juld)
    print('This Julian Day corresponds to {0}/{1}/{2}'.format(gregday[2], gregday[1], gregday[0]))
    return(gregday)


#  ----------------------------------------------------------------------------
def conversion_gregd_juld(day, month, year):
    """Method converting gregorian day into julian day
    
    :rtype: float"""
    #  Petit décalage possible
    julianday = jdcal.gcal2jd(year, month, day)
    juliandayf = julianday[0] + julianday[1] + 0.5
    return juliandayf


if False:
    wmodic = melted.read_dic('wmodic', path_localdata)
    wmodb = melted.read_dic('wmodb', path_localdata)
    argodb = melted.read_dic('argodb', path_localdata)
    #  wmodic = read_wmodic()
    #  wmodb = read_wmstats()
    #  argodb = read_argodb()
