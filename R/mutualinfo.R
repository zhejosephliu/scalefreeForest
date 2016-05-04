# mutual information
mutual.info <- function (x, m, range = NULL) {
  d <- ncol(x)
  K1 <- list()
  K1_val <- matrix(0, nrow = m, ncol = d)
  for (j in c(1:d)) {
    K1[[j]] <- density(x[, j], from = min(x[, j]), to = max(x[, j]), bw = 'nrd', n = m)
    K1_val[, j] <- K1[[j]]$y + 1e-200
  }
  mi <- matrix(0, ncol = d, nrow = d)
  K2 <- list()
  h2 <- rep(0, d)
  for (j in c(1:d)) {
    h2[j] <- bandwidth.nrd2d(x[, j], is1D = F)
    if (h2[j] == 0) {
      stop(paste('Error: not enough distinct values to compute bandwidth in feature', j))
      return(-1)
    }
  }
  h1 <- rep(0, d)
  for (j in c(1:d)) {
    h1[j] <- bandwidth.nrd2d(x[, j], is1D = T)
  }
  for (j in c(1:(d - 1))) {
    for (k in c((j + 1):d)){
      cur_K2 <- kde2d(x[, j], x[, k], n = m, h = c(h2[j], h2[k]) * 4)
      cur_K2_val <- cur_K2$z + 1e-200
      if (is.null(range)) {
        x_range <- max(cur_K2$x) - min(cur_K2$x)
        y_range <- max(cur_K2$y) - min(cur_K2$y)
      } else {
        x_range <- 1
        y_range <- 1
      }
      first_term <- x_range * y_range * sum(cur_K2_val * log(cur_K2_val)) / (m - 1) ^ 2
      second_term <- x_range * y_range * sum(apply(cur_K2_val, 1, sum) / (m - 1) * log(K1_val[, j])) / (m - 1)
      third_term <- x_range * y_range * sum(apply(cur_K2_val, 2, sum) / (m - 1) * log(K1_val[, k])) / (m - 1)
      mi[j, k] <- first_term - second_term - third_term
    }
  }
  return(list(mi = mi))
}
