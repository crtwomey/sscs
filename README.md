
## sscs: Search for Structure in Collective Systems

An implementation of the clustering algorithm introduced in [Twomey et al.
(2018)](https://doi.org/10.1101/362681), provided as an R package. Currently
only implements a hard-clustering variant of this algorithm; the plan is to
include an improved implementation of the soft-clustering variant in the
future. The hard-clustering variant is the most practical for performance
reasons, so it was made the focus of this initial release.

### Compiling

This R package requires Rcpp and RcppEigen (available from CRAN). To compile
from the commandline, cd to the parent directory of the sscs repository and run

```console
$ R CMD build sscs
$ R CMD install sscs_0.1.tar.gz
```

### Usage

The package includes data for a simple example using US state annual average
temperature and precipitation time series data, publicly available from
[NOAA](https://www.ncdc.noaa.gov/cag/statewide/time-series). While simple,
the results are intuitive and easy to visualize.

Univariate example using temperature alone:

```R
# create a new 'sscs' S3 clustering object
sscs <- new_sscs(US_state_temperature)

# run clustering
sscs <- run(sscs, nclusters=12, ncores=16, nreps=400)

# get the cluster assignments
cl <- assignments(sscs)

# convenience function using the 'maps' library for showing the result
plot_us_states_example(cl)
```

Multivariate example combining temperature and precipitation:

```R
# combine measurements
X <- cbind(US_state_temperature, US_state_precipitation)

# js describes how columns of X are grouped into multi-variate variables
nstates <- ncol(US_state_temperature)
js      <- rep(1:nstates, 2)

# run clustering on states
sscs <- new_sscs(X, js)
sscs <- run(sscs, nclusters=12, ncores=16, nreps=400)
cl   <- assignments(sscs)

# show the result
plot_us_states_example(cl)
```

### Credit

Please cite [Twomey et al. (2018)](https://doi.org/10.1101/362681) and include
a link to this repository if you use this code in an academic publication.

Testing, bug reports, and code contributions very welcome.

Copyright (c) 2018, 2019 Colin Twomey.
Shared under a GNU GPLv3 license (see COPYING).

