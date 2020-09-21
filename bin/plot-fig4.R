cli::cli_h1("Plotting Figure 4")
## Setup environment -----
cli::cli_alert_info("Setup environment")

suppressMessages({
  library(tidyverse)
  library(here)
  library(ggbeeswarm)

  source(here("R", "setup.R"))
  source(here("R", "plink.R"))
  source(here("R", "sample.R"))
  source(here("R", "dstat.R"))
  source(here("R", "map.R"))
})

dir.create(here("figs", "fig4"), recursive = TRUE, showWarnings = FALSE)

pluck_pop <- function(dat, ...) {
  filter(dat, ...) %>% pull(pop) %>% unique()
}

## Load data ----
cli::cli_alert_info("Loading data")

plink <- read_plink(here("out", "plink", "DataB"))
fam <- plink$fam %>% left_join(read_sample_info())
pop <- fam %>%
  group_by(pop, group, subgroup, age) %>%
  summarise(across(c(lat, long), mean))

## Compute migration models ------
cli::cli_alert_info("Compute D-statistics models")

native <- pluck_pop(fam, group == "NAT")
amazon <- pluck_pop(fam, subgroup == "Amazon", age == "Contemporan")
wandes <- pluck_pop(fam, subgroup == "West Andes", age == "Contemporan")
south  <- pluck_pop(fam, subgroup == "South", age == "Contemporan")

dstat_models_file <- here("out", "dstat-model.rds")
if (!file.exists(dstat_models_file)) {
  freq <- with(fam, split(sample, pop)) %>%
    discard(~ length(.x) == 0) %>%
    compute_frequency(plink = plink)

  dstat_models <- list(
    model1 = dstat("Mbuti", native, amazon, wandes, freq),
    model2 = dstat("Mbuti", native,  south, wandes, freq),
    model3 = dstat("Mbuti", native, amazon,  south, freq)
  )
  saveRDS(dstat_models, dstat_models_file)
}
dstat_models <- readRDS(dstat_models_file)

## Plots -----
cli::cli_alert_info("Generating plots")

## Fig4. D-stat IDW ------
dstat_idw <- dstat_models %>%
  set_names(c(
    "(A) D(Mbuti,X;Amazon,West Andes)\n(X, Amazon) <---> (X, West Andes)",
    "(B) D(Mbuti,X;South,West Andes)\n(X, South) <---> (X, West Andes)",
    "(C) D(Mbuti,X;Amazon,South)\n(X, Amazon) <---> (X, South)"
  )) %>%
  map(mutate, zv = pmax(pmin(zv, 10), -10)) %>% # Cropping at [-10, 10]
  map(left_join, pop, by = c(x = "pop")) %>%
  idw_america()

dstat_idw_map <- ggplot(pop) +
  geom_raster(data = dstat_idw, aes(x, y, fill = var1.pred)) +
  base_america_map(color = "darkgray", fill = NA) +
  geom_point(aes(long, lat, shape = age)) +
  scale_fill_distiller("Z-scaled D", palette = "RdYlBu") +
  guides(
    shape = guide_legend(
      title.position = "top", direction = "vertical", order = 1
    ),
    color = guide_colorbar(title.position = "top", order = 2),
    fill = guide_colorbar(title.position = "top", order = 2)
  ) +
  facet_grid(~Desc) +
  theme(strip.background = element_blank())

## Save Figure 4
cli::cli_alert_info("Saving Figure 4")

save_plot(
  here("figs", "fig4", "fig4_dstat-idw-map"),
  dstat_idw_map, width = 8, height = 4
)
cli::cli_alert_success("Done!")