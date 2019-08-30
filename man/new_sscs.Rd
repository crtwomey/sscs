\name{new_sscs}
\alias{new_sscs}
\title{ Create a new sscs object. }
\description{
Creates a new sscs S3 object based on data matrix X and variable groupings js.
}
\usage{
  new_sscs(X, js)
}
\arguments{
  \item{X}{ Numeric matrix where rows are observations and columns are variates. }
  \item{js}{ Integer vector with length equal to the number of columns of X, identifying multi-variate groupings of X's columns. }
}
\details{
Currently only permits hard clustering; a soft clustering implementation will be made available in a future release.
}
\value{
Returns a new sscs S3 object to be used with e.g. run(...).
}
\examples{
\dontrun{
# combine measurements
X <- cbind(US_state_temperature, US_state_precipitation)

# describe how columns of X are grouped into
# multi-variate variables.
nstates <- ncol(US_state_temperature)
js      <- rep(1:nstates, 2)

# create a new sscs object to be used with e.g. run(sscs, ...)
sscs <- new_sscs(X, js)

}
}
