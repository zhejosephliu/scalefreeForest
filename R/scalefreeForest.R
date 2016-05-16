#' @import igraph MASS
#' @export

# fit scale-free forest graphical model
scalefreeForest <- function (xtrain, xheld, lambda = seq(0.005, 0.15, 0.005), iter.max = 100, range = NULL, verbose = TRUE) {
    d <- ncol(xtrain)
    if (verbose) {
        cat("Estimating the pairwise mutual information....")
    }
    mi <- mutual.info(xtrain, m = 128, range)$mi
    if (verbose) {
        cat("done.\n")
    }
    wt <- t(mi)[lower.tri(t(mi))]
    g <- graph.full(d)
    E(g)$weight <- wt
    if (verbose) {
        cat("Constructing scale-free forests and computing held-out log-likelihood....\n")
    }
    loglike <- c()
    adj <- list()
    for (lam in c(0, lambda)) {
        forest <- pure.forest(xtrain, xheld, g, wt, bound = 8, prune = F)
        for (i in 1:iter.max) {
            dgr <- degree(forest$seq_g[[d - 1]])
            mdgr <- outer(1 / (dgr + 1e-5), 1 / (dgr + 1e-5), "+")
            wt_cur <- wt - lam * t(mdgr)[lower.tri(t(mdgr))]
            temp <- pure.forest(xtrain, xheld, g, wt_cur, bound = 8, prune = F)
            if (all(get.adjacency(temp$seq_g[[d - 1]]) == get.adjacency(forest$seq_g[[d - 1]]))) {
                if (verbose) {
                    cat(paste('    lambda = ', round(lam, 3), '\t', 'loglike: ', round(max(temp$loglike.fde), 3), '\n', sep = ""))
                }
                loglike <- c(loglike, max(temp$loglike.fde))
                adj[[length(adj) + 1]] <- temp$seq_g[[which.max(temp$loglike.fde)]]
                break
            } else {
                forest <- temp
            }
        }
    }
    if (verbose) {
        cat("done.\n")
    }
    best.loglike <- max(loglike)
    best.adj <- adj[[which.max(loglike)]]
    return(list(loglike = loglike, adj = adj, best.loglike = best.loglike, best.adj = best.adj))
}
