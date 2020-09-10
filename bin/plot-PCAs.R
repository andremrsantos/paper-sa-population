library(tidyverse)
library(here)
library(patchwork)

source(here("R", "setup.R"))
source(here("R", "sample.R"))
source(here("R", "plink.R"))
source(here("R", "pca.R"))

## Setup
theme_set(theme_classic(10, "Helvetica"))

dir.create(here("figs", "sup-fig"), recursive = TRUE, showWarnings = FALSE)
set.seed(05051991)
## Load
dat <- read_plink(here("out", "plink", "DataB"))
dat$fam <- left_join(dat$fam, read_sample_info())

## Amazon PCA
this_pal <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

amz_dat <- filter(dat$fam, study == "This")
amz_pca <- pca(dat, amz_dat$sample)

pca_plot <- amz_dat %>%
  bind_cols(as_tibble(amz_pca$projection)) %>%
  plot_pca_expanded(color = pop) +
  scale_color_manual(values = this_pal)
pve_plot <- plot_pve(amz_pca)

save_plot(
  here("figs", "sup-fig", "sup-fig_amz-pca.pdf"),
  pca_plot / pve_plot +
    plot_layout(heights = c(3,1)) +
    plot_annotation(tag_levels = "A"),
  height = 10, width = 8
)

## Native America PCA
nat_dat <- filter(dat$fam, group == "NAT", age == "Contemporan")
nat_pca <- pca(dat, nat_dat$sample)

pca_plot <- nat_dat %>%
  bind_cols(as_tibble(nat_pca$projection)) %>%
  plot_pca_expanded(color = subgroup) +
  scale_color_manual("Region", values = subgroups_pal)
pve_plot <- plot_pve(nat_pca)

save_plot(
  here("figs", "sup-fig", "sup-fig_nat-pca.pdf"),
  pca_plot / pve_plot +
    plot_layout(heights = c(3,1)) +
    plot_annotation(tag_levels = "A"),
  height = 10, width = 8
)

## Global PCA
global_dat <- filter(dat$fam, age == "Contemporan")
global_pca <- pca(dat, global_dat$sample)

pca_plot <- global_dat %>%
  bind_cols(as_tibble(global_pca$projection)) %>%
  plot_pca_expanded(color = group) +
  scale_color_brewer("Continental Group", palette = "Set2")
pve_plot <- plot_pve(global_pca)

save_plot(
  here("figs", "sup-fig", "sup-fig_global-pca.pdf"),
  pca_plot / pve_plot +
    plot_layout(heights = c(3,1)) +
    plot_annotation(tag_levels = "A"),
  height = 10, width = 8
)

# ## Save pca
# list(
#   amz = amz_dat %>%
#     bind_cols(as_tibble(amz_pca$projection)),
#   nat = nat_dat %>%
#     bind_cols(as_tibble(nat_pca$projection)),
#   global = global_dat %>%
#     bind_cols(as_tibble(global_pca$projection))
#   ) %>%
#   saveRDS(here("out", "pca.rds"))
