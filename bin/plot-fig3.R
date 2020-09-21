cli::cli_h1("Plotting Figure 3")
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
})

dir.create(here("figs", "fig3"), recursive = TRUE, showWarnings = FALSE)

## Load data -----
cli::cli_alert_info("Loading data")

plink <- read_plink(here("out", "plink", "DataB"))
fam <- plink$fam %>% left_join(read_sample_info())
pop <- fam %>%
  group_by(pop, group, subgroup, age) %>%
  summarise(across(c(lat, long), mean))

order <- c(
  "AKW", "ARA", "ARW", "AST", "AWA", "KAY", "ZOE", "WPI", "Surui", "Karitiana",
  "Jabuticabeira", "LapadoSanto", "Laranjal", "Sumidouro"
)

## Computing statistics -----
cli::cli_alert_info("Computing D-statistics for austrolasian signal")

australasian <- unique(filter(fam, group == "OCE", pop != "Hawaiian")$pop)
native       <- unique(filter(fam, group == "NAT", country == "Brazil")$pop)

dstat_australasian_file <- here("out", "dstat-australasian.rds")
if (!file.exists(dstat_australasian_file)) {
  freq <- fam %>%
    filter(pop %in% c(australasian, native, "Mbuti", "Mixe")) %>%
    with(split(sample, pop)) %>%
    discard(~ length(.x) == 0) %>%
    compute_frequency(plink = plink)

  dstat_australasian <- dstat(
    "Mbuti", australasian, "Mixe", native, freq, size = 100
  )
  saveRDS(dstat_australasian, dstat_australasian_file)
}
dstat_australasian <- readRDS(dstat_australasian_file)

## Plotting
cli::cli_alert_info("Generating Figure 3")

ylab <- c(
  "D(Mbuti, Australasin; Mixe, X)\n",
  "(Mixe, Australasian)<--->(X, Australasian)"
)

dstat_scatter <- dstat_australasian %>%
  left_join(pop, by = c(z = "pop")) %>%
  plot_dstat_scatter(z) +
  labs(y = ylab, color = "Age") +
  scale_x_discrete(limits = rev(order)) +
  scale_y_continuous(breaks = c(-10, -5, -2, 0, 2, 5, 10)) +
  coord_flip(ylim = c(-10, 10)) +
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.key.height = unit(1, "lines"),
    legend.key.width = unit(.5, "lines"),
    legend.background = element_blank()
  )

save_plot(
  here("figs", "fig3", "fig3_dstat-australesian"),
  dstat_scatter, width = 4, height = 4.5
)
cli::cli_alert_success("Done!")