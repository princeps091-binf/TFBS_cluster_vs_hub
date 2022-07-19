library(renv)

renv::init()

renv::install("tidyverse")
renv::install("bioc::GenomicRanges")

renv::snapshot()
