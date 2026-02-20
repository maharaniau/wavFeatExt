#' Plot classification results from `classif.wavFeatExt`
#'
#' Create boxplots of cross-validated performance measures (classification
#' error or AUC) for wavelet-based feature sets used in [classif.wavFeatExt()],
#' namely detail and scaling coefficients, and the original/segmented data.
#'
#' @param x The result object returned by [classif.wavFeatExt()], containing at
#'   least the components `CE`, `AUC`, and `method`.
#'
#' @param type Character string specifying which performance measure to plot.
#'   One of `"CE"` for misclassification error, or `"AUC"` for the area under
#'   the ROC curve.
#'
#' @param adjust Logical flag reserved for future extensions of the plotting
#'   function. Currently not used.
#'
#' @param ... Additional graphical parameters passed to [graphics::boxplot()].
#'
#' @details
#' It is assumed that [classif.wavFeatExt()] assigns column names to its `CE`
#' and `AUC` matrices as follows:
#' * Columns beginning with `"D"` (e.g., `"D1"`, `"D2"`, ...) correspond to
#'   wavelet detail coefficients at different scales.
#' * Columns beginning with `"S"` (e.g., `"S1"`, `"S2"`, ...) correspond to
#'   wavelet scaling coefficients at different scales.
#' * The column named `"seg"` corresponds to the original (or segmented) data
#'   matrix.
#'
#' The corresponding columns are extracted and displayed as boxplots, where each
#' box summarises the distribution of cross-validated performance (CE or AUC)
#' across replications for a given feature set.
#'
#' Two horizontal reference lines are added:
#' * A red dashed line at the median performance of the `"seg"` feature set.
#' * A blue dotted line at the best median performance among all feature sets
#'   plotted (maximum median AUC or minimum median CE).
#'
#' @return This function is called for its side effect of producing a plot and
#'   returns `invisible(NULL)`.
#'
#' @author
#' Maharani Ahsani Ummi
#'
#' @seealso
#' [classif.wavFeatExt()],
#' [classif.pcaica()],
#' [wavFeatExt()],
#' [sim.CNA()]
#'
#' @examples
#' ## Generating simulated CNA data
#' set.seed(10)
#' sim.dat10 <- sim.CNA(n.sim = 10)
#'
#' ## Obtain detail and scaling coefficients
#' det.coef.10 <- wavFeatExt(sim.dat10, type = "detail")
#' sca.coef.10 <- wavFeatExt(sim.dat10, type = "scaling")
#'
#' ## Binary response
#' y <- factor(c(rep("Group1", 50), rep("Group2", 50)))
#'
#' ## Classification using Lasso
#' res <- classif.wavFeatExt(sim.dat10, y,
#'                           det.coef.10, sca.coef.10,
#'                           method = "lasso", k = 5)
#'
#' ## Plot classification error for all feature sets
#' plot(res, type = "CE")
#'
#' ## Plot AUC for all feature sets
#' plot(res, type = "AUC")
#'
#' @keywords hplot classification wavelets
#' @export
plot.classif.wavFeatExt <- function(x,
                                    type   = c("CE", "AUC"),
                                    adjust = TRUE,  # reserved
                                    ...) {
  type <- match.arg(type)
  
  ## --- Basic checks ---
  if (!is.list(x) || is.null(x$CE) || is.null(x$AUC)) {
    stop("'x' must be the result object returned by 'classif.wavFeatExt()'.")
  }
  
  CE  <- x$CE
  AUC <- x$AUC
  
  if (is.vector(CE))  CE  <- matrix(CE,  nrow = 1L)
  if (is.vector(AUC)) AUC <- matrix(AUC, nrow = 1L)
  
  if (ncol(CE) != ncol(AUC)) {
    stop("'CE' and 'AUC' in 'x' must have the same number of columns.")
  }
  
  cn <- colnames(AUC)
  if (is.null(cn)) {
    stop("Column names of 'AUC' (and 'CE') must be set in 'classif.wavFeatExt()'.")
  }
  
  ## --- Identifikasi kolom ---
  seg_idx <- which(cn == "seg")
  if (length(seg_idx) != 1L) stop("Could not uniquely identify the 'seg' column in the results.")
  
  # D*, S* (mode per skala)
  d_idx <- grep("^D[0-9]+$", cn)
  s_idx <- grep("^S[0-9]+$", cn)
  
  # ALL (mode all=TRUE). fleksibel: ALL / ALLcoef / all / allcoef
  all_idx <- grep("^ALL($|[A-Za-z0-9_]+$)", cn, ignore.case = TRUE)
  
  ## --- Tentukan kolom yang diplot ---
  if (length(all_idx) >= 1L) {
    # mode all: plot hanya ALL vs seg
    idx <- unique(c(all_idx[1L], seg_idx))
  } else {
    # mode per skala: plot D*, S*, seg (kalau D/S tidak ada, minimal seg tetap)
    idx <- unique(c(d_idx, s_idx, seg_idx))
  }
  
  if (length(idx) < 2L) {
    stop("Not enough columns to plot (need at least 2, e.g., 'ALL' and 'seg' or D*/S* and 'seg').")
  }
  
  ## --- Pilih matriks yang akan diplot ---
  if (type == "AUC") {
    mat  <- AUC[, idx, drop = FALSE]
    main <- paste("AUC:", x$method)
  } else {
    mat  <- CE[, idx, drop = FALSE]
    main <- paste("Classification Error:", x$method)
  }
  
  labels <- colnames(mat)
  
  ## --- Boxplot ---
  boxplot(mat, names = labels, main = main, ...)
  
  ## --- Garis referensi (seg dan best) ---
  seg_col <- which(labels == "seg")
  if (length(seg_col) == 1L) {
    seg_vals  <- mat[, seg_col]
    seg_med   <- stats::median(seg_vals, na.rm = TRUE)
    feat_meds <- apply(mat, 2, stats::median, na.rm = TRUE)
    
    if (type == "AUC") {
      best_val <- max(feat_meds, na.rm = TRUE)   # AUC terbaik = terbesar
    } else {
      best_val <- min(feat_meds, na.rm = TRUE)   # CE terbaik = terkecil
    }
    
    abline(h = seg_med,  lty = 2, col = "red")
    abline(h = best_val, lty = 3, col = "blue")
  }
  
  invisible(NULL)
}
