\name{scalefreeForest}
\alias{scalefreeForest}
\title{Scale-free forest density estimation}
\description{
The main function for scale-free forest density estimation.
}
\usage{
scalefreeForest(xtrain, xheld, lambda = seq(0.005, 0.15, 0.005), iter.max = 100, range = NULL, verbose = TRUE)
}
\arguments{
\item{xtrain}{
An \code{n} by \code{d} data matrix as training set, where \code{n} is the sample size and \code{d} is the dimension.
}
\item{xheld}{
An \code{m} by \code{d} data matrix as held-out set, where \code{m} is the sample size and \code{d} is the dimension.
}
\item{lambda}{
A sequence of positive numbers to control the regularization of the log degree penalty. The default value is \code{seq(0.005, 0.15, 0.005)}.
}
\item{iter.max}{
The maximal number of steps in the iterative reweighted Kruskal's algorithm. The default value is \code{100}.
}
\item{range}{
The range for each of the variables. The default value is \code{NULL}.
}
\item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}
\details{
The training data set is used to construct scale-free forests and the corresponding density estimators, while the held-out data set is then used to determine the optimal scale-free forest by maximizing the held-out log-likelihood.
}
\value{
\item{loglike}{
Maximal held-out log-likelihood for the optimal scale-free forest density estimators corresponding to the regularization parameters.
}
\item{adj}{
Adjacency matrices of the optimal scale-free forests corresponding to the regularization parameters.
}
\item{best.loglike}{
Maximal held-out log-likelihood over all the regularization parameters.
}
\item{best.adj}{
Adjacency matrix of the optimal scale-free forest corresponding to the maximum of the held-out log-likelihood over all the regularization parameters.
}
}
\examples{
library(igraph)

fit <- scalefreeForest(xtrain, xheld)
sf <- fit$best.adj
plot(sf)
}
