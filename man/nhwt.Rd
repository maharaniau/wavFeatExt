\name{nhwt}
\alias{nhwt}
\title{
Non-decimated Haar Wavelet Transform (Legacy)
}
\description{
Performs a non-decimated Haar wavelet transform (stationary wavelet
transform) on one or more one-dimensional signals. This function is
kept for compatibility with older code and is not used in the main
workflow of the package.
}
\usage{
nhwt(data, type = c("detail", "scaling"))
}
\arguments{
  \item{data}{
    Numeric vector, matrix, or data frame containing the data to be
    decomposed. If a matrix or data frame is supplied, rows are
    interpreted as observations and columns as ordered locations.
  }
  \item{type}{
    Character string indicating which coefficients to extract from
    the non-decimated Haar wavelet transform. Use \code{"detail"} for
    detail coefficients or \code{"scaling"} for scaling coefficients.
  }
}
\details{
The function applies a non-decimated (stationary) Haar wavelet
transform to each row of the input. If the number of columns is not a
power of two, the series is extended to the next power of two by
constant padding with value 1 on both sides before computing the
transform. After the transform, the resulting coefficients are
re-aligned and cropped back to the original length.

The underlying transform is computed using the \code{\link[wavethresh]{wd}}
function from the \pkg{wavethresh} package with \code{family = "DaubExPhase"},
\code{filter.number = 1}, and \code{type = "station"} (non-decimated
transform). For each scale, either detail or scaling coefficients can
be extracted via \code{\link[wavethresh]{accessD}} or
\code{\link[wavethresh]{accessC}}, respectively.
}
\value{
A list of length \eqn{J}, where \eqn{J = \lfloor \log_2(n) \rfloor} and
\code{n} is the number of columns in \code{data}. The element at
position \code{j} (for \code{j = 1, \dots, J}) is a numeric matrix of
dimension \code{nrow(data) x n}, containing the detail or scaling
coefficients at scale \code{j} for each observation.
}
\references{
Nason, G. P. (2008).
\emph{Wavelet Methods in Statistics with R}. Springer.
}
\author{
Maharani Ahsani Ummi and Arief Gusnanto
}
\note{
This function is considered legacy and is not used by the main
feature extraction and classification functions in the package.
It is provided for backward compatibility with older analysis
pipelines.
}
\seealso{
\code{\link{wavFeatExt}},
\code{\link{plot.nhwt}},
\code{\link[wavethresh]{wd}},
\code{\link[wavethresh]{accessD}},
\code{\link[wavethresh]{accessC}}
}
\examples{
## Simple example: non-decimated Haar wavelet coefficients
## for a piecewise constant signal
obj <- c(rep(1, 10), rep(3, 20))

## Detail coefficients at all scales
obj.nhwt.det <- nhwt(obj, type = "detail")

length(obj.nhwt.det)      # number of scales
sapply(obj.nhwt.det, dim) # dimensions per scale
}
\keyword{wavelets}
