\name{rescale}
\alias{rescale}
\title{ Rescale an sscs object. }
\description{
Change the quantization scale for an sscs object.
}
\usage{
  rescale(sscs, qscale=NA)
}
\arguments{
  \item{sscs}{ Object of class 'sscs' created by new_sscs(...). }
  \item{qscale}{ Quantization scale (numeric value or NA). }
}
\details{
Changes the quantization scale of the given sscs object to match the numeric value of qscale. If qscale is NA, a sensible default is computed that guarantees the total system redundancy is less than or equal to 1 and each marginal entropy is greater than or equal to 1. Note that computed redundancies for sscs objects are only comparable when they each use the same quantization scale.
}
\value{
Returns the rescaled sscs S3 object.
}
