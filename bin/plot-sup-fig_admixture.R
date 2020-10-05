cli::cli_h1("Plotting Sup. Figure 1/3 - Admixture")
## Setup environment -----
cli::cli_alert_info("Setup environment")

suppressMessages({
  library(tidyverse)
  library(here)

  source(here("R", "setup.R"))
  source(here("R", "sample.R"))
  source(here("R", "plink.R"))
  source(here("R", "admixture.R"))
})

dir.create(here("figs", "sup-fig"), recursive = TRUE, showWarnings = FALSE)

## Load data -----
cli::cli_alert_info("Loading data")

fam <- here("out", "plink", "DataA.fam") %>%
  read_fam() %>%
  left_join(read_sample_info())
dat <- tibble(Ks = 2:10, K = factor(Ks, labels = sprintf("K=%d", Ks))) %>%
  mutate(
    q = here("out", "admixture", "DataA.%d.Q") %>%
      sprintf(Ks) %>%
      map(read_q, fam, subgroup2),
    log = here("out", "admixture", "admixture-K%d.log") %>%
      sprintf(Ks) %>%
      map(read_log)
  ) %>%
  unnest(c(q, log))

## plot full distruct
cli::cli_alert_info("Generating Sup Figure 5 - Admixture")
distruct_plot <-
  plot_distruct(select(dat, -Ks), pop, subgroup2, K = K) +
  scale_fill_manual(values = cbpal) +
  coord_flip() +
  facet_grid(subgroup2~K, scales = "free", space = "free") +
  theme(
    axis.text.y = element_text(size = 5),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 8),
    strip.background.y = element_blank(),
    panel.spacing.x = unit(.4, "lines")
  )
save_plot(
  here("figs", "sup-fig", "sup-fig5_admixture"),
  strip_clip_off(distruct_plot, "r"), width = 8, height = 11
)

cli::cli_alert_info("Generating Sup Figure 6 - Admixture Diagnostics")
full_diagnostic <- dat %>%
  distinct(Ks, cv_error, likelihood) %>%
  rename(`CV Error` = cv_error, Likelihood = likelihood) %>%
  pivot_longer(-Ks) %>%
  ggplot(aes(Ks, value)) +
  geom_line() +
  geom_point() +
  facet_grid(name~., scales = "free", switch = "both") +
  scale_y_continuous("") +
  scale_x_continuous("Puntative Ancestral Components (K)", breaks = 2:10) +
  theme(
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank(),
    strip.placement = "outside"
  )
save_plot(
  here("figs", "sup-fig", "sup-fig6_admixture-diagnostic"),
  full_diagnostic, width = 4, height = 4
)

cli::cli_alert_success("Done!")