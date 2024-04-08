"Run this script to get a list of R packages required, with their current version.

	2024/04/03 @yanisaspic"

suppressPackageStartupMessages({
  library(NCmisc)
})

parse_packages <- function(name) {
  #' Get one or multiple valid package names related to used functions.
  #' 
  #' @param name: a character.
  #' 
  #' @return a vector of characters.
  #' 
  elements <- strsplit(name, split=",")[[1]]
  trim_head <- function(e) {strsplit(e, split="package:")[[1]][2]}
  trim_tail <- function(e) {strsplit(e, split='\\"')[[1]][1]}
  headless_elements <- lapply(X=elements, FUN=trim_head)
  packages <- lapply(X=headless_elements, FUN=trim_tail)
  packages <- unlist(packages)
  return(packages)
}

get_packages.filename <- function(filename) {
  #' Get a vector of packages imported in an R script.
  #' 
  #' @param filename: a character.
  #' 
  #' @return a vector of characters.
  #' 
  source(filename)
  functions <- list.functions.in.file(filename)
  packages <- lapply(X=names(functions), FUN=parse_packages)
  packages <- unlist(packages)
  return(packages)
}

get_packages <- function() {
  #' Get a vector of packages imported in the R scripts of the current directory.
  #' The current directory must be a main directory to prevent issues with imports.
  #' 
  #' @return a vector of characters.
  #' 
  filenames <- list.files(pattern="*\\.R", recursive=TRUE)
  filenames <- unlist(filenames)
  
  packages <- lapply(X=filenames, FUN=get_packages.filename)
  packages <- unlist(packages)
  packages <- packages[!is.na(packages) & !duplicated(packages)]
  packages <- sort(packages)
  return(packages)
}

get_versions <- function(packages) {
  #' Get a list associating packages (names) to their version (values).
  #' 
  #' @param packages: a vector of characters.
  #' 
  #' @return a list associating packages (names) to their version (values).
  #' 
  versions <- lapply(X=packages, FUN=packageVersion)
  names(versions) <- packages
  return(versions)
}

write_versions <- function(versions, path) {
  #' Write a .txt file with the versions of the packages.
  #'
  #' @param versions: a list associating packages (names) to their version (values).
  #' @param path: a character corresponding to the output file.
  #' 
  file_conn <- file(path, "w")
  for (package in names(versions)) {
    writeLines(paste0(package, "==", versions[[package]]), file_conn)
  }
  close(file_conn)
}

report_versions <- function(path) {
  #' Write a .txt file with the versions of the packages imported in the R scripts of the current directory.
  #' The current directory must be a main directory to prevent issues with imports.
  #' 
  #' @param path: a character corresponding to the output file.
  #' 
  packages <- get_packages()
  versions <- get_versions(packages)
  write_versions(versions, path)
}