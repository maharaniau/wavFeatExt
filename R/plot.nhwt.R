#' Plot non-decimated Haar wavelet transform coefficients
#'
#' Plot the non-decimated Haar wavelet transform (NHWT) coefficients
#' produced by [nhwt()] at one or more resolution scales.
#' The function can display either detail or scaling coefficients, and
#' supports global or per-level rescaling for visual comparison.
#'
#' @param x The output from [nhwt()], i.e. a list of matrices
#'   with one element per resolution scale. Each matrix has rows
#'   corresponding to observations and columns to genomic locations
#'   (or data positions). Alternatively, a matrix can be supplied
#'   directly, in which case rows are interpreted as scales and
#'   columns as positions. For multi-observation inputs, only the
#'   first observation is plotted.
#'
#' @param coef Character string indicating the type of coefficients to be
#'   plotted. One of `"detail"` for detail coefficients or
#'   `"scaling"` for scaling coefficients.
#'
#' @param type Character string indicating how the coefficients are rescaled
#'   for plotting when `scale = "all"`. The options are:
#'   \describe{
#'     \item{`"global"`}{
#'       A single scale factor is chosen for the entire plot, based on
#'       the coefficient (across all scales) with the largest absolute
#'       value. This option is useful for comparing coefficients across
#'       different resolution scales.
#'     }
#'     \item{`"by.level"`}{
#'       A separate scale factor is chosen for each resolution scale,
#'       based on the coefficient within that scale that has the largest
#'       absolute value. This option is useful for comparing coefficients
#'       within a given resolution scale.
#'     }
#'   }
#'
#' @param scale Either `"all"` (the default) to plot all resolution scales
#'   in a stacked display, or a single integer specifying the index of
#'   the scale to plot. In the latter case, the chosen scale is shown
#'   on its own.
#'
#' @param ... Further graphical arguments (ignored).
#'
#' @details
#' The function produces a vertical stacked line plot where each
#' resolution scale is shown on a separate horizontal band. Within a
#' scale, each coefficient is drawn as a vertical segment whose height
#' is proportional to its (rescaled) value.
#'
#' If `obj.nhwt` originates from [nhwt()] applied to
#' a data matrix with multiple observations, only the coefficients
#' for the first observation (row) are used. This behaviour is intended
#' for visualising a single time series or copy-number profile.
#'
#' @return This function is called for its side effects of producing a plot.
#'   It returns `invisible(NULL)`.
#'
#' @references
#' Nason, G. P. (2008).
#' *Wavelet Methods in Statistics with R*. Springer.
#'
#' @author
#' Maharani Ahsani Ummi
#'
#' @note
#' This plotting function is primarily intended for exploratory
#' visual inspection of non-decimated Haar wavelet coefficients,
#' for example when assessing multi-scale structure in simulated
#' copy number alteration (CNA) profiles.
#'
#' @seealso
#' [nhwt()],
#' [wavFeatExt()]
#'
#' @examples
#' ## Simple example: single synthetic profile
#' set.seed(1)
#' obj <- c(rep(1, 10), rep(3, 20))
#'
#' ## Non-decimated Haar wavelet transform (detail coefficients)
#' obj.nhwt.det <- nhwt(obj, type = "detail")
#'
#' ## Plot NHWT detail coefficients (all scales, global scaling)
#' plot(obj.nhwt.det, coef = "detail", type = "global")
#'
#' ## Scaling coefficients
#' obj.nhwt.sca <- nhwt(obj, type = "scaling")
#' plot(obj.nhwt.sca, coef = "scaling", type = "global")
#'
#' ## Plot a single scale (scale 1) for detail and scaling coefficients
#' plot(obj.nhwt.det, coef = "detail",  scale = 1)
#' plot(obj.nhwt.sca, coef = "scaling", scale = 1)
#'
#' @keywords hplot wavelets
#' @export
plot.nhwt <- function(x,
                      coef  = c("detail", "scaling"),
                      type  = c("global", "by.level"),
                      scale = "all",
                      ...) {
  obj.nhwt <- x
  
  coef <- match.arg(coef)
  type <- match.arg(type)
  
  ## --- Coerce obj.nhwt into a matrix: rows = scales, cols = locations ---
  if (is.list(obj.nhwt)) {
    n.scales <- length(obj.nhwt)
    if (n.scales < 1L) {
      stop("'obj.nhwt' list is empty.")
    }
    first <- as.matrix(obj.nhwt[[1]])
    n <- ncol(first)
    if (n < 1L) {
      stop("'obj.nhwt' elements must have at least one column.")
    }
    
    coef_mat <- matrix(NA_real_, nrow = n.scales, ncol = n)
    ## For multiple observations, we take the first row of each scale
    for (lev in seq_len(n.scales)) {
      m <- as.matrix(obj.nhwt[[lev]])
      if (ncol(m) != n) {
        stop("All elements of 'obj.nhwt' must have the same number of columns.")
      }
      coef_mat[lev, ] <- as.numeric(m[1, ])
    }
  } else if (is.matrix(obj.nhwt)) {
    coef_mat <- obj.nhwt
    n.scales <- nrow(coef_mat)
    n <- ncol(coef_mat)
    if (n.scales < 1L || n < 1L) {
      stop("'obj.nhwt' matrix must have positive dimensions.")
    }
  } else {
    stop("'obj.nhwt' must be either a list (as returned by nhwt()) or a matrix.")
  }
  
  ## --- Handle 'scale' argument ---
  if (identical(scale, "all")) {
    plot_all_scales <- TRUE
  } else {
    plot_all_scales <- FALSE
    if (!is.numeric(scale) || length(scale) != 1L) {
      stop("'scale' must be \"all\" or a single numeric index.")
    }
    scale <- as.integer(scale)
    if (scale < 1L || scale > n.scales) {
      stop("'scale' must be between 1 and ", n.scales, ".")
    }
  }
  
  ## --- Set up plot labels ---
  main_title <- if (coef == "detail") {
    "Haar NHWT (Detail) Coefficients"
  } else {
    "Haar NHWT (Scaling) Coefficients"
  }
  
  if (plot_all_scales) {
    ## --- Plot all scales ---
    ylim_max <- 2 * n.scales
    plot(1, type = "n",
         xlab = "",
         ylab = "Scale",
         xlim = c(1, n),
         ylim = c(0, ylim_max),
         yaxt = "n",
         main = main_title)
    
    axis(2,
         at     = seq(1, 2 * n.scales - 1, by = 2),
         labels = seq_len(n.scales))
    
    if (type == "global") {
      scale_factor <- max(abs(coef_mat))
      if (scale_factor == 0) {
        scale_factor <- 1
      }
      for (lev in seq_len(n.scales)) {
        vals <- coef_mat[lev, ] / scale_factor
        y0   <- lev * 2 - 1
        for (i in seq_len(n)) {
          lines(c(i, i),
                c(y0, y0 + vals[i]))
        }
      }
    } else if (type == "by.level") {
      for (lev in seq_len(n.scales)) {
        sf <- max(abs(coef_mat[lev, ]))
        if (sf == 0) sf <- 1
        vals <- coef_mat[lev, ] / sf
        y0   <- lev * 2 - 1
        for (i in seq_len(n)) {
          lines(c(i, i),
                c(y0, y0 + vals[i]))
        }
      }
    }
  } else {
    ## --- Plot single scale ---
    plot(1, type = "n",
         xlab = "",
         ylab = "Scale",
         xlim = c(1, n),
         ylim = c(0, 2),
         yaxt = "n",
         main = main_title)
    
    axis(2, at = 1, labels = scale)
    
    sf <- max(abs(coef_mat[scale, ]))
    if (sf == 0) sf <- 1
    vals <- coef_mat[scale, ] / sf
    y0   <- 1  # single scale baseline
    for (i in seq_len(n)) {
      lines(c(i, i),
            c(y0, y0 + vals[i]))
    }
  }
  
  invisible(NULL)
}
