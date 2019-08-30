\name{average_redundancy}
\alias{average_redundancy}
\title{ Compute the average redundancy of an sscs object. }
\description{
Approximates the average redundancy of a cluster solution by computing Gaussian approximations of the normalized total correlation of each cluster.
}
\usage{
  average_redundancy(sscs)
}
\arguments{
  \item{sscs}{ Object of class 'sscs' created by new_sscs(...). }
}
\value{
Returns a numeric scalar.
}
