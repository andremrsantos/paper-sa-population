library(tidyverse)
library(here)
library(Rcpp)

## compile and import read_bed from cpp
Rcpp::sourceCpp(here("src", "read_bed.cpp"))

read_fam <- function(file) {
  cnames <- c("pop", "sample", "father", "mother", "sex", "phenotype")
  readr::read_delim(file, col_names = cnames, col_types = cols(), delim = " ")
}

read_plink <- function(prefix) {
  fam <- read_fam(paste0(prefix, ".fam"))
  bed <- read_bed(paste0(prefix, ".bed"), nrow(fam))
  rownames(bed) <- fam$sample
  list(fam = fam, bed = bed)
}

exec_plink <- function(args) {
  require(sys)
  plink <- c(
    "docker", "run", "--rm", "-w", here(), "-v", paste0(here(), ":", here()),
    "biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1", "plink1.9"
  )
  sys::exec_wait(plink[1], c(plink[-1], args))
}
