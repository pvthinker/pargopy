import numpy as np
import pandas as pd
import argotools as at
import matplotlib.pyplot as plt


zref = at.zref
nz = len(zref)
# let's pretend to have 10 profiles
nprofiles = 10

# let's define random tags
tags = np.random.permutation(range(nprofiles*10))[:nprofiles]


# temp as a numpy array
temp = np.zeros((nprofiles, nz))

# Temp as a DataFrame with tag as index

d = {}
for k, z in enumerate(zref):
    d[z] = temp[:, k]

Temp = pd.DataFrame(d, index=tags)

# accessing all elements of Temps
Temp.loc[:, :]

# accessing all elements at depth 1500
Temp.loc[:, 1500]

# checking that Temp and temp have identical shape
assert np.shape(Temp) == np.shape(temp)
