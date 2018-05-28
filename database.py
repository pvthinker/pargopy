#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 12:39:10 2018

@author: therry

Tools to generate and maintains the processed database

they are high-level routines that rely on smaller modules
"""


header_keys = ['DATA_MODE', 'LONGITUDE', 'LATITUDE', 'JULD', 'FLAG', 'STATUS']


def synchronize_headers():
    """ 
    Propagate the header of each tile onto the global header

    """
    header_global = tools.create_empty_header()
    # pour crÃ©er un empty DataFrame avec les bonnes clefs:
    header_global = pd.DataFrame(columns=header_keys+['TAG'])
    header_global = header_global.set_index('TAG')
    
    for itile in range(ntiles):
        header_tile = tile.read_header(itile)
        header_global =header_global.merge(header_tile,
                                           left_index=True, right_index=True,
                                           how='outer', on=header_keys)
        # warning: because of the overlapping of tiles, some tags
        # appear in several tiles, they need to be exactly identical
        # otherwise merge will complain. In that case, it means that
        # we do someting wrong
    tools.write_header_global(header_global)
   

def read_profile(dac, wmo, iprof=None,
                 header=False, data=False,
                 headerqc=False, dataqc=False,
                 verbose=True):
    """
    :param dac: DAC du profil recherché
    :param wmo: WMO du profil recherché
    :param iprof: Numéro du profil recherché
    :param header: Sélectionne seulement LATITUDE, LONGITUDE et JULD
    :param headerqc: Sélectionne seulement POSITION_QC et JULD_QC
    :param data: Sélectionne TEMP, PSAL et PRES
    :param dataqc: Sélectionne TEMP_QC, PSAL_QC et PRES_QC
    :param verbose: ???
    
    Les valeurs sélectionnée grâce aux arguments passés à la fonction définissent
    la DataFrame que retournera celle-ci.
    
    Basic driver to read the \*_prof.nc data file

    The output is a DataFrame of vectors
    - read one or all profiles read the header (lat, lon, juld) or not
    - read the data or not always return IDAC, WMO, N_PROF, N_LEVELS
    - and DATA_UPDATE (all 5 are int)

    :rtype: DataFrame
    """