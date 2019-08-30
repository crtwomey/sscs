\name{run}
\alias{run}
\title{ Run clustering for sscs object. }
\description{
Run clustering algorithm for a given sscs object. The sscs object will be initialized using the specified number of clusters (nclusters) and any additional parameters (passed in '...', see init(...)). Clustering can be done repeatedly (specified by nreps), returning the best results. Computation will be distributed across the number of compute cores specified by ncores.
}
\usage{
  run(sscs, nclusters=2, nreps=1, ncores=1, max.iterations=100,
      tolerance=0)
}
\arguments{
  \item{sscs}{ Object of class 'sscs' created by new_sscs(...). }
  \item{nclusters}{ Integer number of clusters >= 1. }
  \item{nreps}{ Number of repeated runs (returns the best run). }
  \item{ncores}{ Number of compute cores to distribute repeated runs across. }
  \item{max.iterations}{ Maximum number of assignment update iterations per run. }
  \item{tolerance}{ Numerical precision to use when updating cluster redundancy (see update_assignments) }
  \item{...}{ Additional parameters to pass to init(...). }
}
\details{
Clustering is done by iteratively updating cluster assignments, stopping when assignments are stable (no change between iterations) or max.iterations is reached. Tolerance determines the numerical precision used for calculating cluster redundancies (see update_assignments).
}
\value{
Returns an sscs object with updated cluster assignments, which can be accessed with the assignments(...) function.
}
