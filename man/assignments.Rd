\name{assignments}
\alias{assignments}
\title{ Gets cluster assignments. }
\description{
}
\usage{
  assignment(sscs)
}
\arguments{
  \item{sscs}{ Object of class 'sscs' created by new_sscs(...). }
}
\value{
Returns an integer vector of cluster assignments. For an uninitialized sscs object (one just created by new_sscs(...)) this will simply be a vector of 1's; by default nclusters=1 and all variables belong to that cluster. For an initialized or clustered sscs object (returned by init(sscs, ...) or e.g. run(sscs, ...), respectively) the integers in this vector will vary between 1 and nclusters (sscs$njhat).
}
