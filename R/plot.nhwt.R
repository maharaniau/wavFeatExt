#' Plot NHWT coefficients
#'
#' @param x An object of class \code{nhwt}.
#' @param coef ...
#' @param type ...
#' @param scale ...
#' @param ... Further graphical arguments (ignored).
#'
#' @method plot nhwt
#' @export
plot.nhwt <- function(x,
                      coef  = c("detail", "scaling"),
                      type  = c("global", "by.level"),
                      scale = "all",
                      ...) {
  # biar tidak banyak ubah logika, kita ganti nama lokal:
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
