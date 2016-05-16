# construct forest
pure.forest <- function (xtrain, xheld, g, wt, bound = 8, prune = FALSE) {
  d <- ncol(xtrain)
  n <- nrow(xtrain)
  m <- nrow(xheld)
  h1 <- rep(0, d)
  h2 <- rep(0, d)
  for (j in c(1:d)) {
      h1[j] <- bandwidth.nrd2d(xtrain[, j], is1D = T)
      h2[j] <- bandwidth.nrd2d(xtrain[, j], is1D = F)
  }
  K1 <- matrix(0, nrow = m, ncol = d)
  for (j in 1:d) {
      K1[, j] <- kde1d.new(x = xtrain[, j], h = h1[j], x.new = xheld[, j])$y + 1e-200
  }
  max_edges <- Inf
  seq_g <- list()
  sub_g <- graph.empty(vcount(g))
  sub_g <- as.undirected(sub_g)
  wt_sort_ids <- sort(wt, decreasing = T, index.return = T)$ix
  new_wts <- c()
  loglike.fde <- c()
  cluster_id <- 1:vcount(g)
  cluster_size <- array(1, dim = vcount(g))
  if (max_edges >= ecount(g)) {
    max_edges <- ecount(g)
  }
  if (max_edges == 0) {
    return(sub_g)
  }
  if (ecount(g) > 0) {
    for (i in 1:max_edges) {
      cur_vs <- ends(g, wt_sort_ids[i])
      v1 <- cur_vs[1]
      v2 <- cur_vs[2]
      v1i <- v1
      v2i <- v2
      cur_wt <- wt[wt_sort_ids[i]]
      if (cluster_id[v1i] != cluster_id[v2i] && (!prune || cluster_size[cluster_id[v1i]] + cluster_size[cluster_id[v2i]] <= bound)) {
        cluster_size[cluster_id[v1i]] <- cluster_size[cluster_id[v1i]] + cluster_size[cluster_id[v2i]]
        old_id <- cluster_id[v2i]
        new_id <- cluster_id[v1i]
        for (j in 1:length(cluster_id)) {
          if (cluster_id[j] == old_id) {
            cluster_id[j] <- new_id
          }
        }
        new_wts <- c(new_wts, cur_wt)
        sub_g <- add.edges(sub_g, c(v1, v2))
        seq_g[[length(seq_g) + 1]] <- sub_g
        if(length(loglike.fde) > 0) {
            last.loglike.fde <- loglike.fde[length(loglike.fde)]
        } else {
            last.loglike.fde <- mean(apply(log(K1), 1, sum))
        }
        K2 <- kde2d.new(x = xtrain[, v1], y = xtrain[, v2], h = c(h2[v1], h2[v2]), x.new = xheld[, v1], y.new = xheld[, v2])
        add.loglike.fde <- mean(log(diag(K2$z) + 1e-200) - log(K1[, v1]) - log(K1[, v2]))
        loglike.fde <- c(loglike.fde, last.loglike.fde + add.loglike.fde)
      }
    }
  }
  if (ecount(sub_g) > 0) {
    E(sub_g)$weight <- new_wts
  }
  return(list(seq_g = seq_g, loglike.fde = loglike.fde))
}
