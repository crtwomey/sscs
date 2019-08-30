\name{run_scan}
\alias{run_scan}
\title{ Computes clustering solutions for a range of cluster sizes. }
\description{
Convenience method for scanning a range of cluster sizes. Sequentially executes run(sscs, ...) for each element in the nclusters vector.
}
\usage{
  run_scan(sscs, nclusters=1:10, ...)
}
\arguments{
  \item{sscs}{ Object of class 'sscs' created by new_sscs(...). }
  \item{nclusters}{ Vector of integers, each >= 1. }
  \item{...}{ Additional parameters to pass to run(sscs, ...). }
}
\value{
Returns a list with two named variables, sscs and avg.redundancy.
\describe{
	\item{sscs}{ A list of return values from run(sscs, ...), one for each value in nclusters. }
	\item{avg.redundancy}{ A vector of the average redundancy of each clustering result, calculated using average_redundancy(...). }
}
}
