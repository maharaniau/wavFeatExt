#' Plot classification results from `classif.pcaica`
#'
#' Create boxplots of cross-validated performance measures (classification
#' error or AUC) for the feature sets used in [classif.pcaica()], namely
#' PCA-based, ICA-based, and original/segmented data features.
#'
#' @param x The result object returned by [classif.pcaica()], containing the
#'   components `CE`, `AUC`, and `method`.
#'
#' @param type Character string specifying which performance measure to plot.
#'   One of `"CE"` for misclassification error or `"AUC"` for the area under
#'   the ROC curve.
#'
#' @param d.type Character string indicating which feature types to include
#'   in the plot. The default `"both"` shows results for both PCA-based and
#'   ICA-based features as well as the original/segmented data. Option `"PCA"`
#'   includes only PCA-based features and the original/segmented data, whereas
#'   `"ICA"` includes only ICA-based features and the original/segmented data.
#'
#' @param adjust Logical flag reserved for future extensions of the plotting
#'   function. Currently not used.
#'
#' @param ... Additional graphical parameters passed to [graphics::boxplot()].
#'
#' @details
#' This function assumes that [classif.pcaica()] has assigned column names to
#' its `CE` and `AUC` matrices, with prefixes `"PC"` for PCA-based features,
#' `"I"` for ICA-based features, and the column name `"seg"` for the original
#' (or segmented) data matrix.
#'
#' The corresponding columns are extracted according to `d.type` and displayed
#' as boxplots: each box summarises the distribution of cross-validated
#' performance (CE or AUC) across replications for a given feature set.
#'
#' In addition, two horizontal reference lines are drawn:
#' * The red dashed line shows the median performance of the `"seg"` feature
#'   set (original/segmented data).
#' * The blue dotted line shows the best median performance among all feature
#'   sets plotted (maximum median AUC or minimum median CE).
#'
#' @return This function is used for its side effect of producing a plot and
#'   returns `invisible(NULL)`.
#'
#' @author
#' Maharani Ahsani Ummi
#'
#' @seealso
#' [classif.pcaica()],
#' [classif.wavFeatExt()],
#' [sim.CNA()],
#' [get.pca()],
#' [get.ica()]
#'
#' @examples
#' set.seed(10)
#'
#' ## Simulate CNA data
#' sim.dat10 <- sim.CNA(n.sim = 10)
#'
#' ## Obtain PCA and ICA features
#' pca <- get.pca(sim.dat10, k = 10)
#' ica <- get.ica(sim.dat10, k = 10)
#'
#' ## Binary response
#' y <- factor(c(rep("Group1", 50), rep("Group2", 50)))
#'
#' ## Classification using Lasso
#' res <- classif.pcaica(sim.dat10, y, pca, ica,
#'                       method = "lasso", k = 5)
#'
#' ## Plot classification error (all feature types)
#' plot(res, type = "CE")
#'
#' ## Plot AUC, PCA features + segmented only
#' plot(res, type = "AUC", d.type = "PCA")
#'
#' @keywords hplot classification
#' @export
plot.classif.pcaica <- function(x,
                                type  = c("CE", "AUC"),
                                d.type = c("both", "PCA", "ICA"),
                                adjust = TRUE,  # reserved, not used currently
                                ...) {
  # x: result from classif.pcaica()
  type  <- match.arg(type)
  d.type <- match.arg(d.type)
  
  ## --- Basic checks ---
  if (!is.list(x) || is.null(x$CE) || is.null(x$AUC)) {
    stop("'x' must be the result object returned by 'classif.pcaica()'.")
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
    stop("Column names of 'AUC' (and 'CE') must be set in 'classif.pcaica()'.")
  }
  
  ## --- Identify PCA / ICA / seg columns by name ---
  pc_idx  <- grep("^PC[0-9]+$", cn)
  ic_idx  <- grep("^I[0-9]+$",  cn)
  seg_idx <- which(cn == "seg")
  
  if (length(seg_idx) != 1L) {
    stop("Could not uniquely identify the 'seg' column in the results.")
  }
  
  if (length(pc_idx) == 0L && d.type %in% c("both", "PCA")) {
    warning("No 'PC*' columns found; skipping PCA curves.")
  }
  if (length(ic_idx) == 0L && d.type %in% c("both", "ICA")) {
    warning("No 'I*' columns found; skipping ICA curves.")
  }
  
  idx <- switch(d.type,
                "both" = c(pc_idx, ic_idx, seg_idx),
                "PCA"  = c(pc_idx, seg_idx),
                "ICA"  = c(ic_idx, seg_idx)
  )
  idx <- idx[idx > 0L]  # buang indeks kosong kalau ada
  
  if (length(idx) == 0L) {
    stop("No columns selected to plot for the chosen 'd.type'.")
  }
  
  subset_CE  <- CE[,  idx, drop = FALSE]
  subset_AUC <- AUC[, idx, drop = FALSE]
  labels     <- cn[idx]
  
  ## --- Pilih matriks mana yang dipakai ---
  if (type == "AUC") {
    mat <- subset_AUC
    main_title <- paste("AUC:", x$method)
  } else {
    mat <- subset_CE
    main_title <- paste("Classification Error:", x$method)
  }
  
  ## --- Boxplot ---
  boxplot(mat, names = labels, main = main_title, ...)
  
  ## --- Garis referensi horizontal ---
  # Seg column is always the last of subset (karena selalu kita masukkan terakhir)
  seg_col <- ncol(mat)
  seg_vals <- mat[, seg_col]
  
  seg_med <- stats::median(seg_vals, na.rm = TRUE)
  feat_meds <- matrixStats::colMedians(mat, na.rm = TRUE)
  
  if (type == "AUC") {
    # merah: median AUC seg
    # biru: AUC terbaik (median terbesar)
    best_val <- max(feat_meds, na.rm = TRUE)
    abline(h = seg_med,  lty = 2, col = "red")
    abline(h = best_val, lty = 3, col = "blue")
  } else {
    # merah: median CE seg
    # biru: CE terbaik (median terkecil)
    best_val <- min(feat_meds, na.rm = TRUE)
    abline(h = seg_med,  lty = 2, col = "red")
    abline(h = best_val, lty = 3, col = "blue")
  }
  
  invisible(NULL)
}
