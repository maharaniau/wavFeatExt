\name{wavFeatExt}
\alias{wavFeatExt}
\title{
Feature Extraction Using Non-decimated Haar Wavelet Transform
}
\description{
Performs feature extraction from copy-number (or other one-dimensional)
profiles using the non-decimated (stationary) Haar wavelet transform.
The function computes detail or scaling coefficients at multiple scales
for each observation and organises them into a structured list.
}
\usage{
wavFeatExt(data, type = c("detail", "scaling"))
}
\arguments{
  \item{data}{
    Either a numeric matrix, a list of numeric matrices (with rows
    corresponding to observations and columns to variables or genomic
    locations), or the output from \code{\link{sim.CNA}}. Each matrix
    represents a set of profiles for which wavelet-based features are
    to be extracted.
  }
  \item{type}{
    Character string indicating which type of coefficients to extract
    from the non-decimated wavelet transform. One of \code{"detail"}
    for detail coefficients or \code{"scaling"} for scaling
    coefficients.
  }
}
\details{
For each matrix in \code{data}, the function applies a non-decimated
Haar wavelet transform (via \code{\link{nhwt}}) to each row. This
yields, for each scale, a matrix of wavelet coefficients of the same
dimension as the original data matrix (same number of observations
and locations).

The input \code{data} is first normalised to a list of matrices. If
a single matrix is provided, it is internally wrapped into a list of
length one. When the input is the output of \code{\link{sim.CNA}},
each element of the list corresponds to one simulated data set.

The resulting object is a list of length equal to the number of
matrices in the input. Each element of this list is itself a list of
length equal to the number of available scales. At each scale, a
numeric matrix contains the wavelet coefficients (detail or scaling)
for all observations.
}
\value{
A list of length equal to the number of matrices in \code{data}. For
each data matrix, an element is returned which is a list of length
\eqn{J}, where \eqn{J} is the number of scales determined by the
(non-decimated) wavelet transform. Each sub-list element is a numeric
matrix of the same dimension as the original data matrix, containing
the extracted wavelet coefficients (detail or scaling) for that scale.
}
\references{
Nason, G. P. (2008).
\emph{Wavelet Methods in Statistics with R}. Springer.

Vionanda, et al. (2023).
Prediction of tumour subtypes based on copy number alterations:
Feature extraction for random forest. Submitted.
}
\author{
Maharani Ahsani Ummi and Arief Gusnanto
}
\seealso{
\code{\link{nhwt}},
\code{\link{sim.CNA}},
\code{\link{classif.wavFeatExt}},
\code{\link{classif.pcaica}}
}
\examples{
## Generate simulated CNA data
set.seed(10)
sim.dat <- sim.CNA()   # one simulated data set

## Extract detail coefficients
det.coef <- wavFeatExt(sim.dat, type = "detail")

## Extract scaling coefficients
sca.coef <- wavFeatExt(sim.dat, type = "scaling")

## Not run:
## Generate 10 simulated data sets (this may take some time)
sim.dat10 <- sim.CNA(n.sim = 10)

## Obtain detail coefficients for all simulated data sets
det.coef.10 <- wavFeatExt(sim.dat10, type = "detail")
## End(Not run)
}
\keyword{feature extraction}
\keyword{wavelets}
