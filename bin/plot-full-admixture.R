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
full_distruct <- plot_distruct(select(dat, -Ks), pop, subgroup2, K = K) +
  scale_fill_manual(values = cbpal) +
  facet_grid(K~.)

save_plot(
  here("figs", "sup-fig", "sup-fig_admixture-full"),
  full_distruct, width = 10, height = 8
)

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
  here("figs", "sup-fig", "sup-fig_admixture-diagnostic"),
  full_distruct, width = 10, height = 8
)
