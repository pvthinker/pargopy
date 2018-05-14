# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 13:10:24 2018

@author: herry

Tools used by all the python files to generate the atlas

"""

from __future__ import with_statement
import jdcal
import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pickle as pickle
import param as param

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
                    output['DATA_MODE'] = np.asarray(
                        [c for c in f.variables['DATA_MODE'][idx]])

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


#  ----------------------------------------------------------------------------
def write_dic(name, dic, path_localdata):
    """
    Function used to write each dic used to create .pkl files
    Regroups :
    - write_wmodic, write_wmstats, write_argodb from argodb.py
    - write_argo_filter from research_tools.py
    - write_tile from tile.py

    :rtype: None
    """

    with open('%s/%s.pkl' % (path_localdata, name), 'w') as f:
        pickle.dump(dic, f)

#  ----------------------------------------------------------------------------


def read_dic(name, path_localdata):
    """
    Function used to read each dic used to create .pkl files
    Regroups :
    - read_wmodic, read_wmstats, read_argodb from argodb.py
    - read_argo_filter from research_tools.py
    - read_tile from tile.py

    :rtype: dict"""

    print('read %s.pkl' % name)
    with open('%s/%s.pkl' % (path_localdata, name), 'r') as f:
        dic = pickle.load(f)
    return dic


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
def tile_definition():
    """Define the tiles coordinates, in the form of a vector of lon and
lat + their margins

The tile indexing is

|-----+-----+-----+-----+-----|
| 280 | 281 | 282 | ... | 299 |
|-----+-----+-----+-----+-----|
| ... | ... | ... | ... | ... |
|-----+-----+-----+-----+-----|
|  20 |  21 |  22 | ... |  39 |
|-----+-----+-----+-----+-----|
|   0 |   1 |   2 | ... |  19 |
|-----+-----+-----+-----+-----|


    :rtype: float, float, int, int, float, float

    """

    minlon = -180.
    maxlon = 180.
    # The Arctic is excluded from this procedure
    #
    # TODO: add the Arctic at some point but adopt a more general way
    # of defining tiles (lon x lat rectangles are not suited for the
    # North pole)
    latmin = -80.
    maxlat = 80.

    # number of tiles in longitude and latitude
    # this is hard-coded. We thus have 20x15=300 tiles
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

    # add margins around each tile. This implies some overlapping
    # between tiles. This allows to do statistics at each interior
    # point of the tile, including near the boundaries where margins
    # are large enough to include all nearby profiles.
    marginlat = minmargin / np.cos(latm*np.pi/180)
    marginlon = 2

    return lat, lon, nlat, nlon, marginlat, marginlon


#  ----------------------------------------------------------------------------
def test_tiles(argo, i):
    """Test that tile 'i' is correctly defined with all profiles positions
    within the tile limits, encompassing the margins

    :rtype: None

    """
    # profiles position
    lat = argo['LATITUDE']
    lon = argo['LONGITUDE']

    # tile limits
    latmax = argo['LATMAX_NO_M']
    latmin = argo['LATMIN_NO_M']
    latmargin = argo['MARGINLAT']

    lonmax = argo['LONMAX_NO_M']
    lonmin = argo['LONMIN_NO_M']
    lonmargin = argo['MARGINLON']

    idx1 = np.where(lat > latmax+latmargin)
    idx2 = np.where(lat < latmin-latmargin)
    idx3 = np.where(lon > lonmax+lonmargin)
    idx4 = np.where(lon < lonmin-lonmargin)

    if idx1[0]+idx2[0]+idx3[0]+idx4[0]:
        # the above test is True if the list is not empty
        raise ValueError(
            'There is an error with the dimensions of the tile number %i' % i)


#  ----------------------------------------------------------------------------
def extract_idx_inside_tile(res, argodb):
    """Extract from 'argodb' the list of profiles that are inside the tile
    The tile limits are given by 'res'

       :rtype: dic

    """
    if res['LONMIN_WITH_M'] > res['LONMAX_WITH_M']:
        # special case for tiles extending across the dateline
        idx = np.where((argodb['LATITUDE'] > res['LATMIN_WITH_M'])
                       & (argodb['LATITUDE'] < res['LATMAX_WITH_M'])
                       & ((argodb['LONGITUDE'] > res['LONMIN_WITH_M'])
                          | (argodb['LONGITUDE'] < res['LONMAX_WITH_M'])))
    else:
        # general case
        idx = np.where((argodb['LATITUDE'] > res['LATMIN_WITH_M'])
                       & (argodb['LATITUDE'] < res['LATMAX_WITH_M'])
                       & (argodb['LONGITUDE'] > res['LONMIN_WITH_M'])
                       & (argodb['LONGITUDE'] < res['LONMAX_WITH_M']))

    argo_extract = extract_idx_from_argodb(argodb, idx)
    argo_extract['LATMIN_NO_M'] = res['LATMIN_NO_M']
    argo_extract['LATMAX_NO_M'] = res['LATMAX_NO_M']
    argo_extract['LONMIN_NO_M'] = res['LONMIN_NO_M']
    argo_extract['LONMAX_NO_M'] = res['LONMAX_NO_M']
    argo_extract['MARGINLAT'] = res['MARGINLAT']
    argo_extract['MARGINLON'] = res['MARGINLON']

    return argo_extract


#  ----------------------------------------------------------------------------
def get_idx_from_list_wmo(argodb, wmos):
    """Get the list of profile indices present in argodb that correspond
       to the list of wmos

       :rtype: list of int

    """
    infos = retrieve_infos_from_tag(argodb, argodb['TAG'])
    idx = []
    for w in wmos:
        idx += list(np.where(infos['WMO'] == w)[0])
    return idx


#  ----------------------------------------------------------------------------
def extract_idx_from_argodb(argodb, idx):
    """Return a argodb type dictionnary that is a subset of argodb and
       containing only entries given in idx (list)

       :rtype: dic

    """
    argodb_extract = {}
    for k in argodb.keys():
        argodb_extract[k] = argodb[k][idx]
    return argodb_extract


#  ----------------------------------------------------------------------------
def extract_idx_from_wmostats(wmostats, idx):
    """Return a wmostats type dictionnary that is a subset of wmostats and
       containing only entries given in idx (list)

       :rtype: dic

    """
    wmostats_extract = {}
    keys = wmostats.keys()
    keys.remove('N_WMO')
    for k in keys:
        wmostats_extract[k] = wmostats[k][idx]
    if type(idx) in [int, np.int64]:
        n_wmo = 1
    else:
        n_wmo = len(idx)
    wmostats_extract['N_WMO'] = n_wmo
    return wmostats_extract


#  ----------------------------------------------------------------------------
def conversion_juld_gregd(juld):
    """Method converting julian day into gregorian day

    :rtype: list of int"""
    #  lats, lons, juld = self.reading_variables()
    #  2433282.5000000 corresponds to the Argo origin date
    #  gregday = jdcal.jd2jcal(2433282.500000, juld)
    gregday = jdcal.jd2gcal(2433282.5, juld)
    print('This Julian Day corresponds to {0}/{1}/{2}'.format(
        gregday[2], gregday[1], gregday[0]))

    return(gregday)


#  ----------------------------------------------------------------------------
def conversion_gregd_juld(year, month, day):
    """Method converting gregorian day into julian day

    :rtype: float"""
    #  Petit décalage possible
    julianday = jdcal.gcal2jd(year, month, day)
    juliandayf = julianday[0] + julianday[1]
    return juliandayf - 2433282.5


if False:
    wmodic = read_dic('wmodic', path_localdata)
    wmodb = read_dic('wmodb', path_localdata)
    argodb = read_dic('argodb', path_localdata)
