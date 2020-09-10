library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)

source(here("R", "setup.R"))
source(here("R", "sample.R"))
source(here("R", "plink.R"))
source(here("R", "pca.R"))
source(here("R", "map.R"))

dir.create(here("figs", "fig1"), recursive = TRUE, showWarnings = FALSE)

## Load data ----
plink <- read_plink(here("out", "plink", "DataB"))

fam <- left_join(plink$fam, read_sample_info())
pop <- filter(fam, study == "This") %>%
  group_by(pop) %>%
  summarise(across(c(long, lat), mean))

## Plot population location ----
br_bb <- c(xmin = -58, xmax = -45, ymin = -8, ymax = 4) 
br_sf <- load_br_map()
am_sf <- load_america_map() %>% filter(name != "Brazil")

br_map <- pop %>%
  ggplot() +
  base_america_map(map = am_sf) +
  add_map(map = br_sf) +
  geom_point(aes(long, lat), color = subgroups_pal["Amazon"]) +
  coord_sf(xlim = br_bb[1:2], ylim = br_bb[3:4]) +
  geom_text_repel(aes(long, lat, label = pop), size = 3)

## Plot PCA ---
dat <- filter(fam, study == "This")
pca <- pca(plink, dat$sample)

this_pal <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

pca_plot <- dat %>%
  bind_cols(as_tibble(pca$projection)) %>%
  plot_pca(PC01, PC02, color = pop) +
  scale_color_manual("Population", values = this_pal)

## Mutation Statistics
nat_dat <- here::here("data", "nat_variant_info.tsv.gz") %>%
  readr::read_tsv(col_types = "cicciidccccddd") %>%
  janitor::clean_names() %>%
  dplyr::filter(alt != "*") %>%
  dplyr::mutate(
    ref_af = pmax(kg_af, exac_af, gnmd_af, na.rm = TRUE),
    impact = case_when(
      str_detect(ann_impact, "HIGH") ~ "High",
      str_detect(ann_impact, "MODERATE") ~ "Moderate",
      str_detect(ann_impact, "LOW") ~ "Low",
      str_detect(ann_impact, "MODIFIER") ~ "Modifier",
      TRUE ~ "Modifier"
    ) %>%
      parse_factor(c("High", "Moderate", "Low", "Modifier")),
    known = if_else(is.na(ref_af), "Novel", "Known") %>%
      parse_factor(c("Novel", "Known")),
    freq_class = case_when(
      ac == 1 ~ "Singleton",
      ac < 4 ~ "2-4 Alleles",
      af < 0.1 ~ "10%",
      af < 0.25 ~ "25%",
      TRUE ~ "> 25%"
    ) %>%
      parse_factor(c("Singleton", "2-4 Alleles", "10%", "25%", "> 25%"))
  )

freq_plot <- nat_dat %>%
  ggplot2::ggplot(aes(freq_class, fill = known)) +
  ggplot2::geom_bar() +
  ggplot2::coord_flip() +
  ggplot2::labs(x = "Allele Frequency", y = "Number of Variants", fill = "")
impact_plot <- nat_dat %>%
  ggplot2::ggplot(aes(impact, fill = known)) +
  ggplot2::geom_bar() +
  ggplot2::coord_flip() +
  ggplot2::labs(x = "Impact", y = "Number of Variants", fill = "")

mut_plot <- (impact_plot + guides(fill = "none") + freq_plot) +
  plot_layout() &
  ggplot2::scale_fill_manual(values = c("#fb8072", "#80b1d3"))

## Combine plots
pca_ <- pca_plot +
  geom_point() +
  guides(color = guide_legend(size = 2, ncol = 2)) +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key.height = unit(.75, "line")
  )
mut_ <- mut_plot +
  theme(
    axis.title = element_text(size = 8),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.key.height = unit(.75, "line")
  )

fig1 <- ((br_map + pca_) / mut_) +
  plot_layout(heights = c(3, 1), widths = c(1, 1)) +
  plot_annotation(tag_levels = "A")

save_plot(here("figs", "fig1", "fig-1_amz"), fig1, width = 8, height = 5)