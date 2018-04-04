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
import numpy as np
# import param as param
from netCDF4 import Dataset
import matplotlib.pyplot as plt

daclist = ['aoml', 'bodc', 'coriolis', 'csio',
           'csiro', 'incois', 'jma', 'kma',
           'kordi', 'meds', 'nmdis']

path_argo = '/net/alpha/exports/sciences/data/ARGO/ARGO/201602-ArgoData'

# path where pickle objects are written
path_localdata = '.'

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
    """Retrieve the dac of a wmo"""

    dac = ''
    for d in daclist:
        if wmo in wmodic[d]:
            dac = d
    return dac


def count_wmos(wmodic):
    """Count the total number of wmo in the Argo database base"""

    nwmos = 0
    for dac in daclist:
        nwmos += len(wmodic[dac])
    print('number of wmo in Argo database: %i' % nwmos)
    return nwmos


def count_profiles_in_database(wmostats):
    """Count the total number of profiles in database"""

    nbprofiles = 0
    for nbpr in wmostats['N_PROF']:
        nbprofiles += nbpr
    print('number of profiles in Argo database: %i' % nbprofiles)
    return nbprofiles


def get_profile_file_path(dac, wmo):
    """Return the file path to the *_prof.nc data file"""

    if type(dac) is int:
        dac = daclist[dac]
    filename = '%s/%s/%i/%i_prof.nc' % (path_argo, dac, wmo, wmo)
    return filename


def read_profile(dac, wmo, iprof=None,
                 header=False, data=False,
                 headerqc=False, dataqc=False,
                 verbose=True):
    """Basic driver to read the *_prof.nc data file

    The output is a dictionnary of vectors
    - read one or all profiles read the header (lat, lon, juld) or not
    - read the data or not always return IDAC, WMO, N_PROF, N_LEVELS
    - and DATA_UPDATE (all 5 are int)

    """
    if type(dac) is int:
        dac = daclist[dac]

    filename = get_profile_file_path(dac, wmo)

    if verbose:
        print('/'.join(filename.split('/')[-3:]))

    output = {}

    required_keys = set(['TEMP', 'PSAL', 'PRES'])

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
    """Add the flag to argodb"""

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


def propagate_flag_backward(argodb, subargodb, verbose=True):
    """Update argodb FLAG using subargodb"""

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


def fix_flag_latlonf(argodb):
    """Set a flag error for profiles having bad masked positions

    Bad masked position yield value of 99999. This fix only concerns a
    few profiles from dac='jma'

    """
    idx = np.where(argodb['LONGITUDE'] == 99999.)[0]
    if len(idx) > 0:
        argodb['FLAG'][idx] = 104
    # no return needed because argodb is mutable


def get_tag(kdac, wmo, kprof):
    """Compute the tag number of a profile

    The inverse of get_tag() is retrieve_infos_from_tag()"""

    return (kdac*10000000+wmo)*1000+kprof


def retrieve_infos_from_tag(argodb, tag):
    """Retrieve idac, wmo and iprof from tag (array of int)

    It is the inverse of get_tag()"""

    iprof = tag % 1000
    tag = (tag-iprof) // 1000
    wmo = tag % 10000000
    tag = (tag-wmo) // 10000000
    idac = tag
    output = {'IDAC': idac, 'WMO': wmo, 'IPROF': iprof}
    return output


def plot_location_profiles(argodb):
    """Plot a scatter plot of profiles in argodb

    argodb can be the full database or any subset

    """
    idx = np.where(argodb['FLAG'][:] == 0)[0]
    lon = argodb['LONGITUDE'][idx]
    lat = argodb['LATITUDE'][idx]
    plt.figure()
    plt.plot(lon, lat, '.')
    plt.show()


def plot_wmo_data(dac, wmo):
    """Plot raw 'TEMP' data for dac, wmo"""

    data = read_profile(dac, wmo, data=True)

    varname = 'TEMP'
    plt.figure()
    plt.imshow(data[varname], interpolation='nearest', origin='lower')
    plt.axis('tight')
    plt.title('%s / wmo=%i / %s' % (dac, wmo, varname))
    plt.colorbar()
    plt.show()


def plot_wmos_stats(wmostats):
    """Plot the histogram of number of profiles per number of levels"""

    plt.figure()
    plt.hist(wmostats['N_LEVELS'],
             weights=wmostats['N_PROF'],
             bins=np.arange(0, 2000, 10))
    plt.xlabel('nb of levels')
    plt.ylabel('nb of profiles')
    plt.show()


#  ----------------------------------------------------------------------------
def path_to_wmo():
    """Method creating the first folder containing all the informations of
    the differents profiles known on Argo"""
    path_to_wmo, path_to_argolog = param.data_argotools()
    #  Creating a list with the name of the following repository
    list_folder = []
    list_folder += os.listdir(path_to_wmo)

    #  Sorting the results of list_folder_0 to keep only the good results (without .tar.gz, ...)
    true_dir = [i for i in list_folder if i[-8:] == 'ArgoData']

    for i in range(len(true_dir)):
        list_wmoid_bis = [[], [], [], [], [], [], [], [], [], [], []]
        len_argolog = 2000000
        list_wmoid = np.zeros((len_argolog,), dtype=int)
        size = 0
        list_folder_1 = []
        first_dir_list = []
        dir_file_list = []
        final_dir_list = []
        dac_list = []
        wmo_list = []
        list_folder_1 += os.listdir('%s/%s' % (path_to_wmo, true_dir[i]))
        #  Listing the following repository with the entire adress found
        #  with the first_dir
        current_directory = '%s/%s' % (path_to_wmo, true_dir[i])
        first_dir_list.append(current_directory)

        #  Listing the following datacenter (aoml, bodc, coriolis, ...)

        dac_list += [i for i in list_folder_1 if len(i) <= 8 and i[:] != 'doc']

    for i in range(len(dac_list)):
        current_directory_2 = '%s/%s' % (current_directory, dac_list[i])
        #  list_wmoid[size:size+len(current_directory_2)] = os.listdir('%s' % (current_directory_2))
        #  size += len(current_directory_2)
        list_wmoid_bis[i] += os.listdir('%s' % (current_directory_2))

        for j in range(len(list_wmoid_bis[i])):
            current_directory_3 = "%s/%s" % (current_directory_2, list_wmoid_bis[i][j])
            final_dir_list.append(current_directory_3)
            directory_and_file = "%s/%s/%s" % (current_directory_2, list_wmoid_bis[i][j], list_wmoid_bis[i][j])
            dir_file_list.append(directory_and_file)

    for i, wmo in enumerate(list_wmoid_bis):
        list_wmoid[size:size+len(wmo)] = wmo
        size += len(wmo)

    return (dac_list, wmo_list, final_dir_list, list_wmoid, dir_file_list)

#  ----------------------------------------------------------------------------
    def reading_variables(minlat, maxlat, minlon, maxlon):
        """Salvation of the profiles datas"""
        #  Reading the netCDF file fom ARGO
        if (os.path.isfile(path_to_argolog)):
            f = Dataset(path_to_argolog, "r", format="NETCDF4")
            lats = f.variables['lat'][:]
            lons = f.variables['lon'][:]
            idx = np.where((lats > minlat) & (lats < maxlat) & (lons > minlon) & (lons < maxlon))
            print(len(idx))
            latitudes = f.variables['lat'][idx]
            longitudes = f.variables['lon'][idx]
            prof = f.variables['prof'][idx]
            update = f.variables['update'][idx]
            juld = f.variables['juld'][idx]
            dac = f.variables['dac'][idx]
            wmo = f.variables['wmo'][idx]

            f.close()

        return (latitudes, longitudes, prof, update, juld, dac, wmo)

#  ----------------------------------------------------------------------------
def conversion_juld_gregd(juld):
    """Method converting julian day into gregorian day"""
    #  lats, lons, juld = self.reading_variables()
    #  2433282.5000000 corresponds to the Argo origin date
    gregday = jdcal.jd2jcal(2433282.500000, juld)
    print('This Julian Day corresponds to {0}/{1}/{2}'.format(gregday[2], gregday[1], gregday[0]))
    return(gregday)


#  ----------------------------------------------------------------------------
def conversion_gregd_juld(day, month, year):
    """Method converting gregorian day into julian day"""
    #  Petit décalage possible
    julianday = jdcal.gcal2jd(year, month, day)
    print('This Gregorian Day corresponds to {0}'.format((julianday[0]+julianday[1]+0.5)-2433282.500000))


if False:
    wmodic = read_wmodic()
    wmodb = read_wmstats()
    argodb = read_argodb()
