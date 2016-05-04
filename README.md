# scalefreeForest

This is scalefreeForest R package for scale-free forest density estimation. 

For more information please contact zhe.liu.uchicago@gmail.com

#### Description

We present a nonparametric method for estimating scale-free graphical models. To avoid the usual Gaussian assumption, we restrict the graph to be a forest and build on the work of forest density estimation. The method is motivated from a Bayesian perspective and is equivalent to finding the maximum spanning tree of a weighted graph with a log degree penalty. We solve the optimization problem via a minorize-maximization procedure with Kruskal's algorithm.

#### Installation

To install the [devtools](https://cran.r-project.org/package=devtools) package:

    install.packages("devtools")
    library(devtools)
    install_github("zhejosephliu/scalefreeForest")

#### Usage

Please read the documentation for the main function scalefreeForest.
