from netCDF4 import Dataset
import gr_param as param


def read_profile(dac, wmo, iprof=None):
    """Read one or all (default) profiles of a wmo from the Argo database
    """
    filename = '%s/%s/%i/%i_prof.nc' % (param.argopath,
                                        param.daclist[dac], wmo, wmo)
    print(filename)
    with Dataset(filename, "r", format="NETCDF4") as f:        
        nbprof = len(f.dimensions['N_PROF'])
        if iprof is None:
            idx = range(nbprof)
        else:
            idx = iprof
        temp = f.variables['TEMP'][idx, :]
        temp_qc = f.variables['TEMP_QC'][idx, :]
        sal = f.variables['PSAL'][idx, :]
        sal_qc = f.variables['PSAL_QC'][idx, :]
        pres = f.variables['PRES'][idx, :]
        pres_qc = f.variables['PRES_QC'][idx, :]
        lat = f.variables['LATITUDE'][idx]
        lon = f.variables['LONGITUDE'][idx]
        juld = f.variables['JULD'][idx]

    return temp, sal, pres, temp_qc, sal_qc, pres_qc, lat, lon, juld


# example
temp, sal, pres, temp_qc, sal_qc, pres_qc, lat, lon, juld = read_profile(0, 5900572, iprof=-1)
    
