\name{classif.wavFeatExt}
\alias{classif.wavFeatExt}
\title{
Classification on Wavelet-based Features
}
\description{
Performs classification using several machine learning methods on
features extracted via non-decimated Haar wavelet transformation,
namely detail and scaling coefficients, as well as on the original
(or segmented) data matrix.
}
\usage{
classif.wavFeatExt(data, y, det, sca, method, k, ite = length(data))
}
\arguments{
  \item{data}{
    Either a list of numeric matrices (rows = observations, columns =
    variables), or the output from \code{\link{sim.CNA}}. Each element
    of the list corresponds to one data set or replicate.
  }
  \item{y}{
    Response vector for the classification task, of length equal to
    the number of rows in \code{data[[1]]}. It must be either a binary
    numeric vector coded as \code{0} and \code{1} or a binary factor
    with two levels.
  }
  \item{det}{
    Detail coefficients of \code{data}, typically obtained by
    \code{\link{wavFeatExt}(data, type = "detail")}. It must be a list
    of the same length as \code{data}; for each data set, the
    corresponding element is a list of matrices, one per wavelet scale.
  }
  \item{sca}{
    Scaling coefficients of \code{data}, typically obtained by
    \code{\link{wavFeatExt}(data, type = "scaling")}. It must be a list
    of the same length as \code{data}; for each data set, the
    corresponding element is a list of matrices, one per wavelet scale.
  }
  \item{method}{
    Character string specifying the classification method to use.
    One of \code{"lasso"} for Lasso logistic regression,
    \code{"elnet"} for elastic net, \code{"RF"} for random forest,
    \code{"NN"} for neural network, \code{"PLS"} for partial least
    squares, or \code{"KNN"} for k-nearest neighbours.
  }
  \item{k}{
    Number of folds for k-fold cross-validation. Default is 5.
    The number of observations must be divisible by \code{k}.
  }
  \item{ite}{
    Number of cross-validation replications. The default is
    \code{length(data)}, i.e. one replication per data set in the
    list.
  }
}
\details{
For each replication, the function constructs a k-fold cross-validation
partition of the observations in \code{data[[1]]}. For every feature
set considered (all wavelet detail feature sets, all wavelet scaling
feature sets, and the original/segmented data matrix), the chosen
classification method is trained on the training folds and evaluated
on the held-out fold.

The cross-validated misclassification error (CE) and area under the
ROC curve (AUC) are computed for each fold and then averaged across
folds. This procedure is repeated for \code{ite} replications, and
the resulting averages are stored in matrices for CE and AUC.

The columns of the output matrices correspond to the feature sets
used: the first \code{length(det[[1]])} columns to wavelet detail
feature sets, the next \code{length(sca[[1]])} columns to wavelet
scaling feature sets, and the final column to the original/segmented
data matrix (labelled \code{"seg"}).
}
\value{
A list with the following components:
\item{CE}{
  Numeric matrix of cross-validated misclassification errors.
  Rows correspond to replications (of length \code{ite}) and
  columns to feature sets (wavelet detail features, wavelet
  scaling features, and the original/segmented data).
}
\item{AUC}{
  Numeric matrix of cross-validated areas under the ROC curve.
  The dimensions and column names match those of \code{CE}.
}
\item{method}{
  Character string giving the classification method used.
}
}
\references{
Nason, G. P. (2008).
\emph{Wavelet Methods in Statistics with R}. Springer.

Tibshirani, R. (1996).
Regression shrinkage and selection via the Lasso.
\emph{Journal of the Royal Statistical Society, Series B}
\bold{58}, 267--288.
}
\author{
Maharani Ahsani Ummi and Arief Gusnanto
}
\seealso{
\code{\link{sim.CNA}},
\code{\link{wavFeatExt}},
\code{\link{classif.pcaica}},
\code{\link{nhwt}}
}
\examples{
## Not run:
## Generating simulated CNA data
set.seed(10)
sim.dat2 <- sim.CNA(n.sim = 2)  # two simulated data sets

## Obtain wavelet detail and scaling coefficients
det.coef <- wavFeatExt(sim.dat2, type = "detail")
sca.coef <- wavFeatExt(sim.dat2, type = "scaling")

## Binary response (for example, first 50 vs last 50 samples)
y <- factor(c(rep("Group1", 50), rep("Group2", 50)))

## Perform classification using Lasso
res <- classif.wavFeatExt(sim.dat2, y, det.coef, sca.coef,
                          method = "lasso", k = 5, ite = 2)
## End(Not run)
}
\keyword{classification}
\keyword{wavelets}
