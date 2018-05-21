import numpy as np
import stats
import matplotlib.pyplot as plt
import itertools as it

plt.ion()
date=['2017', '12', '31']
reso=0.5
timeflag='annual'
mode='D'
#itile=228
#tile=stats.date_mode_filter(mode, date, itile)

lon0, lat0 = -40., 34.

lons = lon0+np.arange(-1,1.5,0.5)
lats = lat0+np.arange(-1,1.5,0.5)

nlons=len(lons)
nlats=len(lats)

res=stats.compute_stats_at_zref(mode, date, lons, lats, 0.5)

zref=res['zref']

i0, j0 = 2, 2
# plt.figure(1)
# plt.clf()
# plt.subplot(121)
# plt.plot(res['CT'].T,-zref,'k',alpha=0.1)
# plt.subplot(122)
# plt.plot(res['SA'].T,-zref,'k',alpha=0.1)

v = 'SAmean' #'CTmean'
# plt.plot(res[v][:, j0, i0],-zref)
# plt.plot(res[v][:, j0, i0-2],-zref)
# plt.plot(res[v][:, j0, i0+2],-zref)
# plt.plot(res[v][:, j0-2, i0],-zref)
# plt.plot(res[v][:, j0+2, i0],-zref)

plt.figure(3)
plt.clf()
for k, v in enumerate(['CTbar', 'SAbar', 'CTstd', 'SAstd']):
    plt.subplot(141+k)
    for j, i in it.product(range(nlats),range(nlons)):
        plt.plot(res[v][:, j, i],-zref)
    plt.title(v)

plt.figure(4)
plt.clf()
for k, v in enumerate(['DZstd', 'EAPE']):
    plt.subplot(121+k)
    for j, i in it.product(range(nlats),range(nlons)):
        plt.plot(res[v][:, j, i],-zref)
    plt.title(v)

