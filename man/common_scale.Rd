\name{common_scale}
\alias{common_scale}
\title{ Convert systems to a common scale. }
\description{
Change a list of sscs objects so that they each use the same quantization scale, making their calculated redundancies comparable.
}
\usage{
  common_scale(systems, qscale=NA)
}
\arguments{
  \item{systems}{ List of objects of class 'sscs' created by new_sscs(...). }
  \item{qscale}{ Quantization scale (numeric value or NA). }
}
\details{
Similar to the rescale(...) function, but rescales each sscs object in the systems list to use the same quantization scale. The quantization scale used is either the numeric value provided to qscale or, if qscale is NA, the maximum qscale of the sscs objects in the systems list.
}
\value{
Returns a list of rescaled sscs S3 objects in the same order as they appeared in the systems list.
}
