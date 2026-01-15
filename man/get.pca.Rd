\name{get.pca}
\alias{get.pca}
\title{
Principal Component Analysis (PCA) for a List of Matrices
}
\description{
Applies principal component analysis (PCA) to each matrix in a
collection of data sets and returns cumulative principal components
for use in downstream classification or feature extraction.
}
\usage{
get.pca(x, k, ...)
}
\arguments{
  \item{x}{
    Either a numeric matrix, a list of numeric matrices (with rows
    corresponding to observations and columns to variables), or the
    output from \code{\link{sim.CNA}}. Each matrix will undergo PCA
    separately.
  }
  \item{k}{
    Integer specifying the number of principal components to extract
    from each matrix. Must be at least 2.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[stats]{prcomp}}.
  }
}
\details{
For each matrix in \code{x}, the function calls \code{\link[stats]{prcomp}}
to compute principal component scores. Let \eqn{Z} denote the resulting
score matrix (rows = observations, columns = principal components).

The function then constructs a list of \code{k - 1} matrices per data
set, where the \eqn{j}-th element (for \eqn{j = 2, \dots, k}) contains
the first \eqn{j} principal components, i.e. \code{Z[, 1:j]}. This
cumulative structure mirrors the behaviour of \code{\link{get.ica}} and
provides feature sets with increasing dimensionality (2 components,
3 components, and so on).

When \code{x} is a single matrix or data frame, it is internally
wrapped into a list of length one, and the result is still a list of
lists for consistency.
}
\value{
A list of length equal to the number of matrices in \code{x}. For each
input matrix, an element is returned which is itself a list of length
\code{k - 1}. The \eqn{j}-th element (for \eqn{j = 2, \dots, k}) of
the inner list is a numeric matrix of dimension \code{nrow(x) x j}
containing the first \eqn{j} principal components for that data set.
}
\references{
Jolliffe, I. T. (2002).
\emph{Principal Component Analysis}. 2nd edition. Springer.
}
\author{
Maharani Ahsani Ummi and Arief Gusnanto
}
\note{
The cumulative component representation (2, 3, \dots, \code{k}
components) is designed to be used in combination with
\code{\link{classif.pcaica}} for comparing classification performance
across different feature dimensions.
}
\seealso{
\code{\link{get.ica}},
\code{\link{classif.pcaica}},
\code{\link{sim.CNA}},
\code{\link[stats]{prcomp}}
}
\examples{
## Generate simulated CNA data
set.seed(10)
sim.dat2 <- sim.CNA(n.sim = 2)  # two simulated data sets

## Obtain PCA-based features with 5 components
pca <- get.pca(sim.dat2, k = 5)

length(pca)           # number of data sets
length(pca[[1]])      # number of cumulative PCA feature sets (k - 1)
sapply(pca[[1]], dim) # dimensions for each cumulative feature set
}
\keyword{multivariate}
\keyword{dimension reduction}
