
# coding: utf-8

# In[9]:


import numpy as np
import pandas as pd
import argotools as at
import matplotlib.pyplot as plt

# In[10]:


argo = pd.read_pickle('/home/roullet/data/ARGO/argodb.pkl')


# In[11]:


argo.describe()


# In[12]:

# there are 14 profiles with JULD ~ 2520 in the database they all
# belong to Coriolis.
#
# the minimal value for the others is 17375,
# corresponding to July 28th 1997
datemin = 17375.

# the datemax should never be greater than today
# June 9th 2018 is JULD=24996
datemax = at.today_juld()

argo.loc[np.abs(argo.LONGITUDE) > 180, 'FLAG'] = 31
argo.loc[np.abs(argo.LATITUDE) > 90, 'FLAG'] = 32
argo.loc[argo.JULD > datemax, 'FLAG'] = 33
argo.loc[argo.JULD < datemin, 'FLAG'] = 34

# add columns YEAR, MONTH, DAY, DAYFRAC to argo
res = at.juld2date(argo.JULD)
for k in res.keys():
    argo[k] = res[k]

argo.JULD[argo.FLAG == 0].hist(bins=100)


# In[14]:

# add columns IDAC, WMO, IPROF to argo
res = at.retrieve_infos_from_tag(argo.TAG)
for k in res.keys():
    argo[k] = res[k]


# In[15]:


argo[argo.FLAG == 0].groupby('YEAR').mean()
argo[argo.FLAG == 0].groupby('IDAC').mean()

plt.figure()
for y in range(1997, 2019):
    argo.JULD[(argo.FLAG == 0) & (argo.YEAR == y)].hist(bins=36)
ax = plt.gca()
ax.set_xlim([datemin, datemax])

plt.figure()
argo[argo.JULD > datemax].IDAC.hist()

plt.figure()
argo[argo.JULD < datemin].IDAC.hist()


# number of bad profiles per IDAC
argo[argo.FLAG > 0].groupby('IDAC').FLAG.count()

# number of bad prodiles per FLAG
argo[argo.FLAG > 0].groupby('FLAG').FLAG.count()

# number of bad prodiles per YEAR
argo[argo.FLAG > 0].groupby('YEAR').FLAG.count()

# number of bad prodiles per YEAR and per FLAG
argo[argo.FLAG > 0].groupby(['YEAR', 'FLAG']).FLAG.count()
