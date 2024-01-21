"Run this script to generate artificial scRNA-seq datasets: ZINB-Wave, SPARSim, SymSim.
  See Cao et al. 2021.

	2023/12/20 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(SymSim)
  library(SPARSim)
  library(zimbwave)
})

get_uniform_distribution <- function(n) {
  #' Get a uniform distribution of populations.
  #' 
  #' @param n: number of discrete populations.
  #' 
  #' @return a vector of probability.
  #'
  probs <- rep(1/n, n)
  return(probs)
}

get_geometric_distribution <- function(n, k=1.5) {
  #' Get a geometric distribution of populations.
  #' 
  #' @param n: number of discrete populations.
  #' @param k: factor of multiplication between two populations.
  #' 
  #' @return a vector of probability.
  #' 
  get_prop <- function(x){k**(x)}
  props <- sapply(X=0:(n-1), FUN=get_prop)
  probs <- props / sum(props)
  return(probs)
}

get_random_distribution <- function(n) {
  #' Get a random distribution of populations.
  #' 
  #' @param n: number of discrete populations.
  #' 
  #' @return a vector of probability.
  #' 
  props <- runif(n, min=0.001, max=1/n)
  probs <- props / sum(props)
  return(probs)
}

n_populations <- c("5", "10", "15")
distributions <- c("uniform", "geometric", "random")
proportions <- list()

for (n in n_populations) {
  proportions[[n]] <- list()
  for (dist in distributions) {
    fun <- get(glue("get_{dist}_distribution"))
    prop <- fun(as.numeric(n))
    proportions[[n]][[dist]] <- prop 
  }
}

