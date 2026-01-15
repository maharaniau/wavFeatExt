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
