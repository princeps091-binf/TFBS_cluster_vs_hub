library(renv)

renv::init()

renv::install("tidyverse")
renv::install("furrr")
renv::install("vroom")

renv::install("bioc::GenomicRanges")

renv::snapshot()
