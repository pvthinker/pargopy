Pargopy: Parallel Argo Data Analysis
====================================

**pargopy** is a set of Python modules to perform statistical analyses
 on Argo data using parallel computing.
 
Why pargopy?
------------
 
The `Argo database <https://www.nodc.noaa.gov/argo/>`__ has now
reached a mature stage with almost 2,000,000 profiles over the World
Ocean as of March 2018.

This large number of observations allows to do statistical analyses
everywhere in the ocean within the upper 2,000m. There is enough data
to make 0.5° resolution global atlas.

To make the computation fast **pargopy** offers a set of modules to
split the work into many subtasks that can be computed using several
cores on a cluster.


Which statistics?
------------------

Several layers of statistics:

- Eulerian means of temperature, salinity, (in-situ) density
- standard deviations, including isopycnal displacement and EAPE
- higher-order statistics (skewness)
- covariance
- etc


Licence
-------

Not yet defined
