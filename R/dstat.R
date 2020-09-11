library(tidyverse)
library(Rcpp)
library(matrixStats)
library(ggbeeswarm)

Rcpp::sourceCpp(here("src", "d.cpp"))

as_matrix <- function(lst)
  dplyr::bind_cols(lst) %>% as.matrix()

row_sum   <- function(cols, data)
  matrixStats::rowSums2(data, cols = cols, na.rm = TRUE)

allele_count <- function(indv, mtx) {
  matrixStats::colSums2(mtx, rows = indv, na.rm = TRUE)
}

allele_total <- function(indv, mtx) {
  total_na <- matrixStats::colCounts(mtx, rows = indv, value = NA_integer_)
  2 * (length(indv) - total_na)
}

compute_count <- function(plink, pops) {
  count <- pops %>% 
    map(match, rownames(plink$bed)) %>%
    map(allele_count, mtx = plink$bed) %>%
    as_matrix()
  colnames(count) <- names(pops)
  return(count)
}

compute_total <- function(plink, pops) {
  total <- pops %>% 
    map(match, rownames(plink$bed)) %>%
    map(allele_total, mtx = plink$bed) %>%
    as_matrix()
  colnames(total) <- names(pops)
  return(total)
}

compute_frequency <- function(plink, pops) {
  allele_count <- compute_count(plink, pops)
  allele_total <- compute_total(plink, pops)
  
  return(allele_count / allele_total)
}

## Compute D-statistic -----------
as_idx <- function(names, w, x, y, z = NULL) {
  idx <- cbind(match(w, names), match(x, names), match(y, names))
  if (!is.null(z)) { idx <- cbind(idx, match(z, names)) }
  return(idx)
}

dstat <- function(w, x, y, z, freq, size = 100, fn = jackknife_ds) {
  list(w=w, x=x, y=y, z=z) %>%
    cross_df() %>%
    filter(w != x, w != y, w != z, x != y, x != y, y != z) %>%
    mutate(
      jk = fn(freq, as_idx(colnames(freq), w, x, y, z), size),
      d = map_dbl(jk, "mean"),
      zv = map_dbl(jk, "z"),
      n = map_dbl(jk, "n"),
      nb = map_dbl(jk, "n_block")
    )
}

f3_stat <- function(w, x, y, ac, an, size = 100, fn = jackknife_f3s) {
  list(w = w, x = x, y = y) %>%
    cross_df() %>%
    filter(w != x, w != y, x != y) %>%
    mutate(
      jk = fn(ac, an, as_idx(colnames(ac), w, x, y), size),
      d = map_dbl(jk, "mean"),
      zv = map_dbl(jk, "z"),
      n = map_dbl(jk, "n"),
      nb = map_dbl(jk, "n_block")
    )
}

plot_dstat_scatter <- function(dat, x = x, color = age) {
  ggplot(dat, aes(reorder({{x}}, {{color}}), zv)) +
    geom_hline(yintercept = c(-3, 0, 3), linetype = "dashed", color = "firebrick") +
    geom_quasirandom(size = 1, shape = 1) +
    stat_summary(aes(color = {{color}}), fun.data = median_hilow) +
    coord_flip(ylim = c(-15, 15)) +
    labs(x = "Population") +
    scale_color_brewer(palette = "Set1") +
    theme(
      legend.position ="bottom",
      strip.text.y = element_text(angle = 0, hjust = 0),
      strip.background = element_blank()
    )
}
