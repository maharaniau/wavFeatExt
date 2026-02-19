\name{plot.classif.wavFeatExt}
\alias{plot.classif.wavFeatExt}
\title{
Plot classification results from \code{classif.wavFeatExt}
}
\description{
Create boxplots of cross-validated performance measures (classification
error or AUC) for wavelet-based feature sets used in
\code{\link{classif.wavFeatExt}}, namely detail and scaling
coefficients, and the original/segmented data.
}
\usage{
plot.classif.wavFeatExt(x,
                        type   = c("CE", "AUC"),
                        adjust = TRUE, ...)
}
\arguments{
  \item{x}{
    The result object returned by \code{\link{classif.wavFeatExt}},
    containing at least the components \code{CE}, \code{AUC}, and
    \code{method}.
  }
  \item{type}{
    Character string specifying which performance measure to plot.
    One of \code{"CE"} for misclassification error, or \code{"AUC"}
    for the area under the ROC curve.
  }
  \item{adjust}{
    Logical flag reserved for future extensions of the plotting
    function. Currently not used.
  }
  \item{\dots}{
    Additional graphical parameters passed to \code{\link{boxplot}}.
  }
}
\details{
It is assumed that \code{\link{classif.wavFeatExt}} assigns column
names to its \code{CE} and \code{AUC} matrices as follows:
\itemize{
  \item Columns beginning with \code{"D"} (e.g., \code{"D1"}, \code{"D2"}, \dots)
        correspond to wavelet detail coefficients at different scales.
  \item Columns beginning with \code{"S"} (e.g., \code{"S1"}, \code{"S2"}, \dots)
        correspond to wavelet scaling coefficients at different scales.
  \item The column named \code{"seg"} corresponds to the original
        (or segmented) data matrix.
}

The corresponding columns are extracted and displayed as boxplots,
where each box summarises the distribution of cross-validated
performance (CE or AUC) across replications for a given feature set.

Two horizontal reference lines are added:
\itemize{
  \item A red dashed line at the median performance of the
        \code{"seg"} feature set.
  \item A blue dotted line at the best median performance among all
        feature sets plotted (maximum median AUC or minimum median CE).
}
}
\value{
This function is called for its side effect of producing a plot and
returns \code{invisible(NULL)}.
}
\author{
Maharani Ahsani Ummi and Arief Gusnanto
}
\seealso{
\code{\link{classif.wavFeatExt}},
\code{\link{classif.pcaica}},
\code{\link{wavFeatExt}},
\code{\link{sim.CNA}}
}
\examples{
## Not run:
## Generating simulated CNA data
set.seed(10)
sim.dat10 <- sim.CNA(n.sim = 10)

## Obtain detail and scaling coefficients
det.coef.10 <- wavFeatExt(sim.dat10, type = "detail")
sca.coef.10 <- wavFeatExt(sim.dat10, type = "scaling")

## Binary response
y <- factor(c(rep("Group1", 50), rep("Group2", 50)))

## Classification using Lasso
res <- classif.wavFeatExt(sim.dat10, y,
                          det.coef.10, sca.coef.10,
                          method = "lasso", k = 5)

## Plot classification error for all feature sets
plot(res, type = "CE")

## Plot AUC for all feature sets
plot(res, type = "AUC")
## End(Not run)
}
\keyword{hplot}
\keyword{classification}
\keyword{wavelets}
