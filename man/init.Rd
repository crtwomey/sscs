\name{init}
\alias{init}
\title{ Initialize an sscs object. }
\description{
Specify the number of clusters and initialize cluster membership assignments for a given sscs object.
}
\usage{
  init(sscs, nclusters=2, balanced=FALSE)
}
\arguments{
  \item{sscs}{ Object of class 'sscs' created by new_sscs(...). }
  \item{nclusters}{ Integer number of clusters >= 1. }
  \item{balanced}{ Boolean determining whether or not to bias the initial cluster memberships to be evenly distributed across clusters.  Set to TRUE only if performance is an issue, otherwise use the default (FALSE) to keep the initial memberships unbiased. }
}
\value{
Returns an initialized sscs S3 object that can be iteratively updated using update_assignment(...).
}
