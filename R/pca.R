# devtools::install_github("gabraham/flashpca/flashpcaR")
library(flashpcaR)
library(matrixStats)
library(tidyverse)
library(ggforce)

pca <- function(plink, samples, pcs = 20, freq_cutoff = 0.01) {
  rows <- match(samples, rownames(plink$bed))
  cols <- which(colMeans2(plink$bed, rows = rows, na.rm = TRUE)/2 > freq_cutoff)

  scaled <- flashpcaR::scale2(plink$bed[rows, cols])
  pca <- flashpcaR::flashpca(scaled, ndim = pcs, stand = "none")
  colnames(pca$projection) <- sprintf("PC%02d", 1:pcs)
  return(pca)
}

plot_pca <- function(dat, x, y, ..., shape = 1) {
  ggplot2::ggplot(dat, aes({{x}}, {{y}}, ...)) +
    ggplot2::geom_point(shape = shape, alpha = .75) +
    ggplot2::theme(
      legend.justification = c(1, 0),
      legend.position = c(0.99, 0.01),
      panel.grid = element_blank(),
      legend.key.height=unit(0.7,"line")
    )
}

plot_pca_expanded <- function(dat, range = PC01:PC04, color = pop) {
  ggplot2::ggplot(dat, aes(.panel_x, .panel_y, color = {{color}})) +
    ggplot2::geom_point() +
    ggforce::facet_matrix(vars({{range}}), switch = 'both') +
    ggplot2::theme(
      panel.border = element_rect(fill = NA),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside"
    )
}

plot_pve <- function(pca) {
  pve <- pca$values / sum(pca$values)
  ggplot2::qplot(x = seq_along(pve), y = pve, geom = "col") +
    ggplot2::scale_x_continuous(expand = expansion(add = .5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "PCs", y = "Percentage of Variance Explained")
}