#! /usr/bin/env Rscript --vanilla

## Install requirements
installed_packages <- rownames(installed.packages())
required <- c(
  "tidyverse", "here", "furrr", "patchwork", "ggforce", "ggrepel", "ggbeeswarm",
  "tidygraph", "ggraph", "remotes", "Rcpp", "RcppArmadillo", "matrixStats",
  "sf", "rnaturalearth", "sp", "gstat", "raster", "sys", "ape", "cli"
)
install.packages(required[!required %in% installed_packages])
## flashpcaR
if (!"flashpcaR" %in% installed_packages)
  remotes::install_github("gabraham/flashpca/flashpcaR")
## admixtools
if (!"admixtools" %in% installed_packages)
  remotes::install_github("uqrmaie1/admixtools")
## ggtree
if (!"ggtree" %in% installed_packages)
  BiocManager::install("ggtree")
