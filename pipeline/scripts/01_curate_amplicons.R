#!/usr/bin/env Rscript
# library(here)    # optional, for paths
source("R/functions1.R")
source("R/helpers.R")

args <- commandArgs(trailingOnly=TRUE)

# Snakemake automatically provides 'snakemake' object when using script:
input_file  <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# Example: run analysis
data <- read.csv(input_file)
res <- analyze_data(data)
saveRDS(res, output_file)
