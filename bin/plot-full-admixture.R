library(tidyverse)
library(here)

source(here("R", "setup.R"))
source(here("R", "sample.R"))
source(here("R", "plink.R"))
source(here("R", "admixture.R"))

strip_clip_off <- function(plot, strip = c("r", "l", "t", "b")) {
  strip <- paste0("strip-", match.arg(strip))
  q <- ggplotGrob(plot)
  for (i in which(grepl(strip, q$layout$name)))
    q$grobs[[i]]$layout$clip <- "off"
  patchwork::wrap_elements(q)
}

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
distruct <- plot_distruct(select(dat, -Ks), pop, subgroup2, K = K) +
  scale_fill_manual(values = cbpal) +
  facet_grid(K~subgroup2, scales = "free", space = "free")

full_distruct <- distruct +
  facet_grid(K~subgroup2, scales = "free", space = "free") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
    strip.text.x = element_text(angle = 90, hjust = 0, size = 8),
    strip.background.x = element_blank()
  )
save_plot(
  here("figs", "sup-fig", "sup-fig_admixture-full"),
  strip_clip_off(full_distruct, "t"), width = 10, height = 8
)

vert_distruct <- distruct +
  facet_grid(subgroup2~K, scales = "free", space = "free") +
  coord_flip() +
  theme(
    axis.text.y = element_text(size = 5),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 8),
    strip.background.y = element_blank(),
    panel.spacing.x = unit(.4, "lines")
  )
save_plot(
  here("figs", "sup-fig", "sup-fig_admixture-vertical"),
  strip_clip_off(vert_distruct, "r"), width = 8, height = 11
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
  full_diagnostic, width = 4, height = 4
)
