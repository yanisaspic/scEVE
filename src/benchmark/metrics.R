"Functions used to compute metrics for the benchmark.

	2024/01/21 @yanisaspic"

suppressPackageStartupMessages({
  library(aricode)
})

get_ARI <- function(query, truth) {
  #' Compute the ARI score from two vectors.
  #' 
  #' @param query: a vector of labels being evaluated.
  #' @param truth: a vector of ground truths.
  #' 
  #' @return a numeric.
  #' 
  result <- ARI(query, truth)
  return(result)
}

get_NMI <- function(query, truth) {
  #' Compute the NMI score from two vectors.
  #' 
  #' @param query: a vector of labels being evaluated.
  #' @param truth: a vector of ground truths.
  #' 
  #' @return a numeric.
  #' 
  result <- NMI(query, truth)
  return(result)
}

get_time <- function(query, truth) {
  #' Get the computation time.
  #' 
  #' @param query: a vector of identical times (in seconds).
  #' @param truth: filler for SummarizedBenchmark.
  #' 
  #' @return a numeric (in seconds).
  #' 
  result <- max(query)
  return(result)
}

get_memory <- function(query, truth) {
  #' Get the peak RAM.
  #' 
  #' @param query: a vector of identical peak RAMs (in Mb).
  #' @param truth: filler for SummarizedBenchmark.
  #' 
  #' @return a numeric (in Mb).
  #' 
  result <- max(query)
  return(result)
}