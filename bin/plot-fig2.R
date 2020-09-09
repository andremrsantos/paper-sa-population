library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)
library(ggforce)

source(here("R", "setup.R"))
source(here("R", "sample.R"))
source(here("R", "plink.R"))
source(here("R", "pca.R"))
source(here("R", "map.R"))
source(here("R", "admixture.R"))

dir.create(here("figs", "fig2"), recursive = TRUE, showWarnings = FALSE)

## Load data ----
plink <- read_plink(here("out", "plink", "DataB"))

fam <- left_join(plink$fam, read_sample_info())
pop <- filter(fam, group == "NAT") %>%
  group_by(pop, subgroup, age) %>%
  summarise(across(c(long, lat), mean))

qfm <- read_fam(here("out", "plink", "DataA.fam")) %>% left_join(read_sample_info())
qnm <- sprintf("K=%d", 5:7)
qdt <- here("out", "admixture", "DataA.%d.Q") %>%
  sprintf(5:7) %>%
  set_names(qnm) %>%
  map_dfr(read_q, qfm, subgroup2, .id = "K") %>%
  filter(sample %in% fam$sample, group != "OCE") %>%
  mutate(K = parse_factor(K, qnm))

## Plot samples map
america_map <- ggplot(pop) +
  base_america_map() +
  geom_point(aes(long, lat, shape = age, color = subgroup)) +
  geom_text_repel(aes(long, lat, label = pop), size = 2) +
  scale_color_manual("Region", values = subgroups_pal) +
  guides(
    color = guide_legend(title.position = "top", direction = "vertical"),
    shape = guide_legend(title.position = "top", direction = "vertical")
  )

## Plot distruct
america_distruct <-
  plot_distruct(qdt, subgroup2, subgroup2, pop, K = K) +
  scale_fill_manual(values = cbpal) +
  facet_grid(K~., scales = "free", space = "free")

## Plot scatterpie
america_scatterpie <- qdt %>%
  filter(K == "K=7", !is.na(long), !is.na(lat)) %>%
  pivot_longer(
    starts_with("K") & where(is.numeric),
    names_to = "cluster",
    values_to = "rate"
  ) %>%
  group_by(pop, age, cluster) %>%
  summarise(across(c(long, lat, rate), mean)) %>%
  ggplot() +
  base_america_map() +
  geom_arc_bar(
    aes(x0 = long, y0 = lat, r0 = 0, r = 2, amount = rate, fill = cluster),
    stat = "pie", color = NA, show.legend = FALSE
  ) +
  scale_fill_manual(values = cbpal) +
  facet_grid(age~.) +
  theme(strip.background = element_blank())

## Plot PCA
nat_dat <- fam %>% filter(group == "NAT", age == "Contemporan")
nat_pca <- pca(plink, nat_dat$sample)

nat_pca_plot <- bind_cols(nat_dat, as_tibble(nat_pca$projection)) %>%
  plot_pca(PC01, PC02, color = subgroup) +
  geom_point() +
  scale_color_manual("Region", values = subgroups_pal, guide = "none")

## Combine plots
fig2 <- (
  (america_map + america_scatterpie + plot_layout(widths = c(2, 1))) /
  (nat_pca_plot + america_distruct + plot_layout(widths = c(1, 3)))
) +
  plot_layout(heights = c(4, 1)) +
  plot_annotation(tag_levels = "A")
ggsave(
  here("figs", "fig2", "fig-2_america.pdf"),
  fig2, width = 8, height = 9, useDingbats = FALSE
)
