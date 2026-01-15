\name{wavFeatExt-package}
\alias{wavFeatExt-package}
\alias{wavFeatExt}
\docType{package}
\title{wavFeatExt: Wavelet-based Feature Extraction for Copy-number Alteration Data}

\description{
Utility functions for non-decimated Haar wavelet transformation of copy-number alteration (CNA) profiles, feature extraction, and supervised learning workflows.
}

\details{
The package provides tools for simulating CNA profiles, performing segmentation, extracting wavelet-derived features (detail and scaling coefficients) at multiple scales, and fitting classifiers using wavelet features or PCA/ICA features.
}

\author{
\packageAuthor{wavFeatExt}

Maintainer: \packageMaintainer{wavFeatExt}
}

\seealso{
\itemize{
  \item \code{\link{sim.CNA}} for simulating CNA profiles.
  \item \code{\link{wavFeatExt}} for wavelet-based feature extraction.
  \item \code{\link{classif.wavFeatExt}} for classification using wavelet features.
}
}
\keyword{package}
