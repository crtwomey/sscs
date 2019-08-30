\name{sscs-package}
\alias{sscs-package}
\alias{sscs}
\docType{package}
\title{
Search for structure in collective systems via redundancy partitioning
}
\description{
An implementation of the clustering algorithm introduced in \href{https://doi.org/10.1101/362681}{Twomey et al. (2018)}, provided as an R package. Currently only an implementation of the hard-clustering variant is provided; the plan is to include an improved implementation of the soft-clustering variant in a future release. The hard-clustering variant is by far the most practical for performance reasons, so it was made the focus of this initial release.
}
\details{
	\tabular{ll}{
		Package: \tab sscs\cr
		Type: \tab Package\cr
		Version: \tab 0.1\cr
		Date: \tab 2019-08-28\cr
		License: \tab GPL-3\cr
	}
}
\author{
	Colin Twomey

	Maintainer: Colin Twomey <crtwomey@gmail.com>
}
\references{
Twomey et al. (2018) Searching for structure in collective systems. \url{https://doi.org/10.1101/362681}
}
\examples{
#
# A simple example using US state average temperature time
# series data, publicly available from NOAA here:
#     https://www.ncdc.noaa.gov/cag/statewide/time-series
#
# The results are intuitive and easy to visualize with a map.
#
\dontrun{
# univariate example
# ------------------

# create a new 'sscs' S3 clustering object
sscs <- new_sscs(US_state_temperature)

# run clustering
sscs <- run(sscs, nclusters=12, ncores=16, nreps=400)

# get the cluster assignments
cl <- assignments(sscs)

# convenience function for showing the result
plot_us_states_example(cl)


# multivariate example
# --------------------

# combine measurements
X <- cbind(US_state_temperature, US_state_precipitation)

# js describes how columns of X are grouped into
# multi-variate variables.
nstates <- ncol(US_state_temperature)
js      <- rep(1:nstates, 2)

# run clustering on states
sscs <- new_sscs(X, js)
sscs <- run(sscs, nclusters=12, ncores=16, nreps=400)
cl   <- assignments(sscs)

# show the result
plot_us_states_example(cl)

}
}
