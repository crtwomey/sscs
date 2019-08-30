\name{update_assignment}
\alias{update_assignment}
\title{ Update cluster assignment. }
\description{
Updates cluster assignment for variable istar, where istar is one of the (possibly multi-variate) variables being clustered. Generally it is easiest to use the run(sscs, ...) interface for generating clusters, but this method can be used for more granular control over the iterative updating of the cluster solution.
}
\usage{
  update_assignment(sscs, istar, tolerance=0)
}
\arguments{
  \item{sscs}{ Object of class 'sscs' created by new_sscs(...). }
  \item{istar}{ Integer index of variable to update (must be <= sscs$nj). }
  \item{tolerance}{ 0, 1, or 2 corresponding to high, medium, and low precision estimates of cluster redundancy. }
}
\details{
Low precision may be useful for speeding up calculations when the number of elements, nj, is very large. Otherwise it can actually be slower than just using the higher precision method.
}
\value{
Returns an updated sscs object.
}
