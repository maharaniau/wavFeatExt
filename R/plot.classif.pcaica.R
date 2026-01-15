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
