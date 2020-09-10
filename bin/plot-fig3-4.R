library(tidyverse)
library(here)
library(ggbeeswarm)

source(here("R", "setup.R"))
source(here("R", "plink.R"))
source(here("R", "sample.R"))
source(here("R", "dstat.R"))
source(here("R", "map.R"))

## Load data ----
plink <- read_plink(here("out", "plink", "DataB"))

dat <- left_join(plink$fam, read_sample_info())
fam <- plink$fam %>% left_join(read_sample_info())

pop_sample <- fam %>%
  with(c(split(sample, pop), split(sample, subgroup))) %>%
  discard(~ length(.x) == 0)
pop_freq   <- compute_frequency(plink, pop_sample)

pop_order <- c(
  "AKW", "ARA", "ARW", "AST", "AWA", "KAY", "PTJ", "WPI", "Surui", "Karitiana",
  "Jabuticabeira", "LapadoSanto", "Laranjal", "Sumidouro"
)

## Compute migration models ------
## Listing populations per region
native <- unique(filter(dat, group == "NAT")$pop)
amazon <- unique(filter(dat, subgroup == "Amazon", age == "Contemporan")$pop)
wandes <- unique(filter(dat, subgroup == "West Andes", age == "Contemporan")$pop)
south  <- unique(filter(dat, subgroup == "South", age == "Contemporan")$pop)

dstat_models_file <- here("out", "dstat-model.rds")
if (!file.exists(dstat_models_file)) {
  dstat_models <- list(
    model1 = dstat("Mbuti", native, amazon, wandes, pop_freq),
    model2 = dstat("Mbuti", native,  south, wandes, pop_freq),
    model3 = dstat("Mbuti", native, amazon,  south, pop_freq)
  )
  saveRDS(dstat_models, dstat_models_file)
} else {
  dstat_models <- readRDS(dstat_models_file)
}

## Plots ------
## SupFig. D-stat Amazon v West Andes -----
dstat_info <- dat %>%
  group_by(pop, group, subgroup, age) %>%
  summarise(across(c(lat, long), mean))

dstat_scatter <- dstat_models[["model1"]] %>%
  left_join(dstat_info, by = c(z = "pop")) %>%
  filter(y == "Mixe", x != z, z %in% pop_order) %>%
  plot_dstat_scatter() +
  facet_grid(subgroup ~ ., scales = "free_y", space = "free") +
  labs(
    y = "D(Mbuti, X; Amazon, WestAndes)\n(X, Amazon) <---> (X, WestAndes)",
    color = "Age"
  )

# dir.create(here("figs", "sup-fig"), recursive = TRUE, showWarnings = FALSE)
# save_plot(
#   here("figs", "sup-fig", "sup-fig_dstat-amazon-v-wandes"),
#   dstat_scatter, width = 4, height = 5
# )

## SupFig. D-stat migration models -----
dstat_scatter <- bind_rows(
  mutate(dstat_models[["model1"]],  Y = "Amazon", Z = "West Andes"),
  mutate(dstat_models[["model2"]], Y = "South", Z = "West Andes"),
  mutate(dstat_models[["model3"]], Y = "Amazon", Z = "South")
  ) %>%
  mutate(D = sprintf("D(Mbuti, X; %s, %s)\n(X,%s)<--->(X,%s)", Y, Z, Y, Z)) %>%
  left_join(dstat_info, by = c(x = "pop")) %>%
  plot_dstat_scatter() + 
  facet_grid(D~subgroup, scales = "free", space = "free") +
  coord_cartesian(ylim = c(-15, 15)) +
  theme(
    strip.text.x = element_text(angle = 90, hjust = 0),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# dir.create(here("figs", "sup-fig"), recursive = TRUE, showWarnings = FALSE)
# save_plot(
#   here("figs", "sup-fig", "sup-fig_dstat-scatter"),
#   dstat_scatter, width = 8, height = 4
# )

## Fig4. D-stat IDW ------
dstat_idw <- dstat_models %>%
  set_names(c(
    "D(Mbuti,X;Amazon,West Andes)\n(X, Amazon) <---> (X, West Andes)",
    "D(Mbuti,X;South,West Andes)\n(X, South) <---> (X, West Andes)",
    "D(Mbuti,X;Amazon,South)\n(X, Amazon) <---> (X, South)"
  )) %>%
  map(mutate, zv = pmax(pmin(zv, 10), -10)) %>% # Cropping at [-10, 10]
  map(left_join, dstat_info, by = c(x = "pop")) %>%
  idw_america()

dstat_idw_map <- ggplot(dstat_info) +
  geom_raster(data = dstat_idw, aes(x, y, fill = var1.pred)) +
  base_america_map(color = "darkgray", fill = NA) +
  geom_point(aes(long, lat, shape = age)) +
  scale_fill_distiller("Z-scaled D", palette = "RdYlBu") +
  guides(
    shape = guide_legend(title.position="top", direction = "vertical", order = 1),
    color = guide_colorbar(title.position="top", order = 2),
    fill = guide_colorbar(title.position="top", order = 2)
  ) +
  facet_grid(~Desc) +
  theme(strip.background = element_blank())

dir.create(here("figs", "fig4"), recursive = TRUE, showWarnings = FALSE)
save_plot(
  here("figs", "fig4", "fig4_dstat-idw-map"),
  dstat_idw_map, width = 8, height = 4
)

## Fig3. D-stat Australasian -----
australasian <- unique(filter(dat, group == "OCE", pop != "Hawaiian")$pop)
native       <- unique(filter(dat, group == "NAT", country == "Brazil")$pop)
native_anc   <- unique(filter(dat, group == "NAT", age == "Ancient")$pop)

dstat_australasian_file <- here("out", "dstat-australasian.rds")
if (!file.exists(dstat_australasian_file)) {
  dstat_australasian <- dstat(
    "Mbuti", australasian, c("Mixe", native_anc), c(native, australasian),
    pop_freq, size = 100
  )
  saveRDS(dstat_australasian, dstat_australasian_file)
} else {
  dstat_australasian <- readRDS(dstat_australasian_file)
}

dstat_scatter <- dstat_australasian %>%
  left_join(dstat_info, by = c(z = "pop")) %>%
  filter(!x %in% c("Hawaiian"), y == "Mixe", x != z, z %in% pop_order) %>%
  plot_dstat_scatter(z) +
  labs(
    y = "D(Mbuti, Australasin; Mixe, X)\n(Mixe, Australasian)<--->(X, Australasian)",
    color = "Age"
  )  +
  scale_x_discrete(limits = rev(pop_order)) +
  scale_y_continuous(breaks = c(-10, -5, -2, 0, 2, 5, 10)) +
  coord_flip(ylim = c(-10, 10)) +
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.key.height = unit(1, "lines"),
    legend.key.width = unit(.5, "lines"),
    legend.background = element_blank()
  )

dir.create(here("figs", "fig3"), recursive = TRUE, showWarnings = FALSE)
save_plot(
  here("figs", "fig3", "fig3_dstat-australesian"),
  dstat_scatter, width = 4, height = 4.5
)