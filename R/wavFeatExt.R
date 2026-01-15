wavFeatExt <- function(data, type = c("detail", "scaling")) {
  # Feature extraction using non-decimated Haar wavelet transform (NHWT)
  # data: matrix, list of matrices, or output from sim.CNA()
  # type: "detail" or "scaling"
  
  type <- match.arg(type)
  
  ## Normalise input to a list of matrices
  if (is.matrix(data)) {
    data_list <- list(data)
  } else if (is.data.frame(data)) {
    data_list <- list(as.matrix(data))
  } else if (is.list(data)) {
    data_list <- lapply(data, function(x) {
      if (!is.matrix(x)) as.matrix(x) else x
    })
  } else {
    stop("'data' must be a matrix, data frame, or a list of matrices.")
  }
  
  n.sim <- length(data_list)
  if (n.sim == 0L) {
    stop("'data' must contain at least one data matrix.")
  }
  
  ## Apply NHWT to each data matrix
  res <- vector("list", n.sim)
  for (j in seq_len(n.sim)) {
    res[[j]] <- nhwt(data_list[[j]], type = type)
  }
  
  res
}
