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

get_AMI <- function(query, truth) {
  #' Compute the AMI score from two vectors.
  #' 
  #' @param query: a vector of labels being evaluated.
  #' @param truth: a vector of ground truths.
  #' 
  #' @return a numeric.
  #' 
  result <- AMI(query, truth)
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
  result <- log10(as.numeric(query[1], units="secs")+1)
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
  result <- log10(query[1]+1)
  return(result)
}

frame_scEVE_metrics <- function(scEVE_results, ground_truth) {
  #' Frame the results of a scEVE analysis to the metrics table.
  #' 
  #' @param scEVE_results: a list with three items: 'labels', 'peakRAM' and 'memory'.
  #' @param ground_truth: a vector of expected labels.
  #' 
  #' @return a data.frame with three columns: 'label', 'value' and 'performanceMetric'.
  #' 
  metrics <- c("AMI", "ARI", "NMI")
  for (met in metrics) {
    scEVE_results[[met]] <- get(glue("get_{met}"))(scEVE_results$labels, ground_truth)
  }
  scEVE_results$labels <- NULL
  result <- stack(scEVE_results)
  colnames(result) <- c("value", "performanceMetric")
  return(result)
}

add_scEVE_metrics <- function(metrics, scEVE_results, ground_truth) {
  #' Add the results of the scEVE analysis to the metrics table.
  #' 
  #' @param metrics: a data.frame with three columns: 'label', 'value' and 'performanceMetric'.
  #' @param scEVE_results: a list with three items: 'labels', 'peakRAM' and 'memory'.
  #' @param ground_truth: a vector of expected labels.
  #' 
  #' @return a data.frame with three columns: 'label', 'value' and 'performanceMetric'.
  #' 
  get_resolution <- function(label){nchar(label) - 1}
  resolutions <- sapply(X=scEVE_results$labels, FUN=get_resolution)
  for (j in 1:max(resolutions)) {
    data <- scEVE_results
    data$labels <- sapply(X=scEVE_results$labels, FUN=substr, start=1, stop=1+j)
    
    scEVE_metrics <- frame_scEVE_metrics(scEVE_results=data, ground_truth=ground_truth)
    scEVE_metrics$label <- glue("{j}.scEVE")
    scEVE_metrics$highlight <- "yes"
    
    metrics <- rbind(metrics, scEVE_metrics)
  }
  return(metrics)
}
