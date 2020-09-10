library(tidyverse)
library(here)

source(here("R", "setup.R"))
source(here("R", "sample.R"))
source(here("R", "plink.R"))
source(here("R", "admixture.R"))

## Setup
dir.create(here("figs", "sup-fig"), recursive = TRUE, showWarnings = FALSE)

## Load data
fam <- here("out", "plink", "DataA.fam") %>%
  read_fam() %>%
  left_join(read_sample_info())

nms <- sprintf("K=%d", 2:10)
dat <- here("out", "admixture", "DataA.%d.Q") %>%
  sprintf(2:10) %>%
  set_names(nms) %>%
  map_dfr(read_q, fam, subgroup2, .id = "K") %>%
  mutate(K = parse_factor(K, nms))

## plot full distruct
full_distruct <- plot_distruct(dat, pop, subgroup2, K = K) +
  scale_fill_manual(values = cbpal) +
  facet_grid(K~.)

save_plot(
  here("figs", "sup-fig", "sup-fig_admixture-full"),
  full_distruct, width = 10, height = 8
)
