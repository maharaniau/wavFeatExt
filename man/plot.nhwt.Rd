\name{plot.nhwt}
\alias{plot.nhwt}
\title{
Plot non-decimated Haar wavelet transform coefficients
}
\description{
Plot the non-decimated Haar wavelet transform (NHWT) coefficients
produced by \code{\link{nhwt}} at one or more resolution scales.
The function can display either detail or scaling coefficients, and
supports global or per-level rescaling for visual comparison.
}
\usage{
plot.nhwt(obj.nhwt,
          coef  = c("detail", "scaling"),
          type  = c("global", "by.level"),
          scale = "all")
}
\arguments{
  \item{obj.nhwt}{
    The output from \code{\link{nhwt}}, i.e. a list of matrices
    with one element per resolution scale. Each matrix has rows
    corresponding to observations and columns to genomic locations
    (or data positions). Alternatively, a matrix can be supplied
    directly, in which case rows are interpreted as scales and
    columns as positions. For multi-observation inputs, only the
    first observation is plotted.
  }
  \item{coef}{
    Character string indicating the type of coefficients to be
    plotted. One of \code{"detail"} for detail coefficients or
    \code{"scaling"} for scaling coefficients.
  }
  \item{type}{
    Character string indicating how the coefficients are rescaled
    for plotting when \code{scale = "all"}. The options are:
    \describe{
      \item{\code{"global"}}{
        A single scale factor is chosen for the entire plot, based on
        the coefficient (across all scales) with the largest absolute
        value. This option is useful for comparing coefficients across
        different resolution scales.
      }
      \item{\code{"by.level"}}{
        A separate scale factor is chosen for each resolution scale,
        based on the coefficient within that scale that has the largest
        absolute value. This option is useful for comparing coefficients
        within a given resolution scale.
      }
    }
  }
  \item{scale}{
    Either \code{"all"} (the default) to plot all resolution scales
    in a stacked display, or a single integer specifying the index of
    the scale to plot. In the latter case, the chosen scale is shown
    on its own.
  }
}
\details{
The function produces a vertical stacked line plot where each
resolution scale is shown on a separate horizontal band. Within a
scale, each coefficient is drawn as a vertical segment whose height
is proportional to its (rescaled) value.

If \code{obj.nhwt} originates from \code{\link{nhwt}} applied to
a data matrix with multiple observations, only the coefficients
for the first observation (row) are used. This behaviour is intended
for visualising a single time series or copy-number profile.
}
\value{
This function is called for its side effects of producing a plot.
It returns \code{invisible(NULL)}.
}
\references{
Nason, G. P. (2008).
\emph{Wavelet Methods in Statistics with R}. Springer.
}
\author{
Maharani Ahsani Ummi and Arief Gusnanto
}
\note{
This plotting function is primarily intended for exploratory
visual inspection of non-decimated Haar wavelet coefficients,
for example when assessing multi-scale structure in simulated
copy number alteration (CNA) profiles.
}
\seealso{
\code{\link{nhwt}},
\code{\link{wavFeatExt}}
}
\examples{
## Simple example: single synthetic profile
set.seed(1)
obj <- c(rep(1, 10), rep(3, 20))

## Non-decimated Haar wavelet transform (detail coefficients)
obj.nhwt.det <- nhwt(obj, type = "detail")

## Plot NHWT detail coefficients (all scales, global scaling)
plot.nhwt(obj.nhwt.det, coef = "detail", type = "global")

## Scaling coefficients
obj.nhwt.sca <- nhwt(obj, type = "scaling")
plot.nhwt(obj.nhwt.sca, coef = "scaling", type = "global")

## Plot a single scale (scale 1) for detail and scaling coefficients
plot.nhwt(obj.nhwt.det, coef = "detail",  scale = 1)
plot.nhwt(obj.nhwt.sca, coef = "scaling", scale = 1)
}
\keyword{hplot}
\keyword{wavelets}
