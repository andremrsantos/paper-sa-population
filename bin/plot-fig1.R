cli::cli_h1("Plotting Figure 1")
## Setup environment -----
cli::cli_alert_info("Setup environment")

suppressMessages({
  library(tidyverse)
  library(here)
  library(patchwork)
  library(janitor)

  source(here("R", "setup.R"))
  source(here("R", "sample.R"))
  source(here("R", "plink.R"))
  source(here("R", "pca.R"))
  source(here("R", "map.R"))
})

dir.create(here("figs", "fig1"), recursive = TRUE, showWarnings = FALSE)

this_pal <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

## Load data ----
cli::cli_alert_info("Loading data")
plink <- read_plink(here("out", "plink", "DataB"))

fam <- left_join(plink$fam, read_sample_info())
pop <- filter(fam, study == "This") %>%
  group_by(pop) %>%
  summarise(across(c(long, lat), mean), .groups = "drop")

## Generate plots ------
cli::cli_alert_info("Preparing panels")

## Plot population location ----
plot_box   <- c(xmin = -58, xmax = -45, ymin = -8, ymax = 4)
brazil_sf  <- load_br_map()
america_sf <- load_america_map() %>% filter(name != "Brazil")
rivers_sf  <- ne_download(
  type = "rivers_lake_centerlines", scale = 50, category = "physical",
  returnclass = "sf"
  ) %>%
  ## Crop to plotting region
  st_crop(plot_box) %>%
  ## Collapse rivers
  group_by(label) %>%
  summarise()

plot_map <- ggplot(pop) +
  ## Plot maps
  base_america_map(map = america_sf) +
  add_map(map = brazil_sf) +
  geom_sf(data = rivers_sf, color = "steelblue", alpha = .5) +
  ## Add river labels
  geom_sf_text(
    aes(label = label, geometry = geometry),
    data = filter(rivers_sf, label != "Xingu"), size = 2.5,
    color = "steelblue", nudge_y = .8, nudge_x = -.4,
  ) +
  geom_sf_text(
    aes(label = label, geometry = geometry),
    data = filter(rivers_sf, label == "Xingu"), size = 2.5,
    color = "steelblue", nudge_y = -.75, nudge_x = -.25, hjust = 1,
  ) +
  ## Indicate populations location
  geom_point(aes(long, lat, color = pop)) +
  geom_text(
    aes(long, lat, label = pop), size = 3,
    nudge_x = .1, nudge_y = .25, hjust = 0
  ) +
  coord_sf(xlim = plot_box[1:2], ylim = plot_box[3:4], expand = FALSE) +
  scale_color_manual(values = this_pal)

## Plot PCA -----
dat <- filter(fam, study == "This")
pca <- pca(plink, dat$sample)

plot_pca <- dat %>%
  bind_cols(as_tibble(pca$projection)) %>%
  plot_pca(PC01, PC02, color = pop) +
  scale_color_manual("Population", values = this_pal) +
  geom_point() +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key.height = unit(.75, "line")
  )

## Plot Mutation Statistics -----
## Load variants info
nat_dat <- here("data", "nat_variants.tsv.gz") %>%
  read_tsv(col_types = "cicccdccc") %>%
  mutate(
    known  = parse_factor(is_known, c("Novel", "Known")),
    impact = parse_factor(
      allele_impact, c("High", "Moderate", "Low", "Modifier")
    ),
    freq_class = parse_factor(
      frequency_class, c("Singleton", "2-4 Alleles", "10%", "25%", "> 25%")
    )
  )

plot_by_freq <- ggplot(nat_dat, aes(freq_class, fill = known)) +
  geom_bar() +
  coord_flip() +
  labs(x = "Allele Frequency", y = "Number of Variants", fill = "")
plot_by_impact <- ggplot(nat_dat, aes(impact, fill = known)) +
  geom_bar(show.legend = FALSE) +
  coord_flip() +
  labs(x = "Impact", y = "Number of Variants", fill = "")

plot_snv <- (plot_by_impact + plot_by_freq) &
  scale_fill_manual(values = c("#fb8072", "#80b1d3")) &
  theme(
    axis.title = element_text(size = 8),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.key.height = unit(.75, "line")
  )

## Combine plots
cli::cli_alert_info("Saving Figure 1")

map_ <- wrap_elements(plot = plot_map + theme(plot.margin = margin()))
save_plot(
  here("figs", "fig1", "fig-1_amz"),
  ((map_ + plot_pca) / plot_snv) +
    plot_layout(heights = c(4, 1)) +
    plot_annotation(tag_levels = "A"),
  width = 8, height = 5
)

cli::cli_alert_success("Done!")