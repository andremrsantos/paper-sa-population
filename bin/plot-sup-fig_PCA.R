cli::cli_h1("Plotting Sup. Figure 4-6 - PCA")
cli::cli_alert_info("Setup environment")

suppressMessages({
  library(tidyverse)
  library(here)
  library(patchwork)

  source(here("R", "setup.R"))
  source(here("R", "sample.R"))
  source(here("R", "plink.R"))
  source(here("R", "pca.R"))
})

dir.create(here("figs", "sup-fig"), recursive = TRUE, showWarnings = FALSE)

## Load
cli::cli_alert_info("Load data")
dat <- read_plink(here("out", "plink", "DataB"))
fam <- left_join(dat$fam, read_sample_info())

## Amazon PCA
cli::cli_alert_info("Generating Sup Figure 1 - Amazon")
this_pal <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

amz_dat <- filter(fam, study == "Present study")
amz_pca <- pca(dat, amz_dat$sample)

pca_plot <- amz_dat %>%
  bind_cols(as_tibble(amz_pca$projection)) %>%
  plot_pca_expanded(color = pop) +
  scale_color_manual("Population", values = this_pal)
pve_plot <- plot_pve(amz_pca)

save_plot(
  here("figs", "sup-fig", "sup-fig1_amz-pca"),
  pca_plot / pve_plot +
    plot_layout(heights = c(3, 1)) +
    plot_annotation(tag_levels = "A"),
  height = 10, width = 8
)

## Native America PCA
cli::cli_alert_info("Generating Sup Figure 4 - Contemporan Native American")
nat_dat <- filter(fam, group == "NAT", age == "Contemporan")
nat_pca <- pca(dat, nat_dat$sample)

pca_plot <- nat_dat %>%
  bind_cols(as_tibble(nat_pca$projection)) %>%
  plot_pca_expanded(color = subgroup) +
  scale_color_manual("Region", values = subgroups_pal)
pve_plot <- plot_pve(nat_pca)

save_plot(
  here("figs", "sup-fig", "sup-fig4_nat-pca"),
  pca_plot / pve_plot +
    plot_layout(heights = c(3, 1)) +
    plot_annotation(tag_levels = "A"),
  height = 10, width = 8
)

## Global PCA
cli::cli_alert_info("Generating Sup Figure 2 - Whole world")
global_dat <- filter(fam, age == "Contemporan")
global_pca <- pca(dat, global_dat$sample)

pca_plot <- global_dat %>%
  bind_cols(as_tibble(global_pca$projection)) %>%
  plot_pca_expanded(color = group) +
  scale_color_brewer("Continental Group", palette = "Set2")
pve_plot <- plot_pve(global_pca)

save_plot(
  here("figs", "sup-fig", "sup-fig2_global-pca"),
  pca_plot / pve_plot +
    plot_layout(heights = c(3, 1)) +
    plot_annotation(tag_levels = "A"),
  height = 10, width = 8
)

cli::cli_alert_success("Done!")