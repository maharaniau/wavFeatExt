\name{plot.classif.pcaica}
\alias{plot.classif.pcaica}
\title{
Plot classification results from \code{classif.pcaica}
}
\description{
Create boxplots of cross-validated performance measures (classification
error or AUC) for the feature sets used in
\code{\link{classif.pcaica}}, namely PCA-based, ICA-based, and
original/segmented data features.
}
\usage{
plot.classif.pcaica(x,
                    type  = c("CE", "AUC"),
                    d.type = c("both", "PCA", "ICA"),
                    adjust = TRUE, ...)
}
\arguments{
  \item{x}{
    The result object returned by \code{\link{classif.pcaica}},
    containing the components \code{CE}, \code{AUC}, and \code{method}.
  }
  \item{type}{
    Character string specifying which performance measure to plot.
    One of \code{"CE"} for misclassification error or \code{"AUC"}
    for the area under the ROC curve.
  }
  \item{d.type}{
    Character string indicating which feature types to include in the
    plot. The default \code{"both"} shows results for both PCA-based
    and ICA-based features as well as the original/segmented data.
    Option \code{"PCA"} includes only PCA-based features and the
    original/segmented data, whereas \code{"ICA"} includes only
    ICA-based features and the original/segmented data.
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
This function assumes that \code{\link{classif.pcaica}} has assigned
column names to its \code{CE} and \code{AUC} matrices, with prefixes
\code{"PC"} for PCA-based features, \code{"I"} for ICA-based features,
and the column name \code{"seg"} for the original (or segmented) data
matrix.

The corresponding columns are extracted according to \code{d.type}
and displayed as boxplots: each box summarises the distribution of
cross-validated performance (CE or AUC) across replications for a
given feature set.

In addition, two horizontal reference lines are drawn:
\itemize{
  \item The red dashed line shows the median performance of the
        \code{"seg"} feature set (original/segmented data).
  \item The blue dotted line shows the best median performance among
        all feature sets plotted (maximum median AUC or minimum
        median CE).
}
}
\value{
This function is used for its side effect of producing a plot and
returns \code{invisible(NULL)}.
}
\author{
Maharani Ahsani Ummi and Arief Gusnanto
}
\seealso{
\code{\link{classif.pcaica}},
\code{\link{classif.wavFeatExt}},
\code{\link{sim.CNA}},
\code{\link{get.pca}},
\code{\link{get.ica}}
}
\examples{
## Not run:
set.seed(10)
## Simulate CNA data
sim.dat10 <- sim.CNA(n.sim = 10)

## Obtain PCA and ICA features
pca <- get.pca(sim.dat10, k = 10)
ica <- get.ica(sim.dat10, k = 10)

## Binary response
y <- factor(c(rep("Group1", 50), rep("Group2", 50)))

## Classification using Lasso
res <- classif.pcaica(sim.dat10, y, pca, ica,
                      method = "lasso", k = 5)

## Plot classification error (all feature types)
plot(res, type = "CE")

## Plot AUC, PCA features + segmented only
plot(res, type = "AUC", d.type = "PCA")
## End(Not run)
}
\keyword{hplot}
\keyword{classification}
