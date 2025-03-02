#! Rscript --vanilla

cli::cli_h1("Plotting Sup. Figure 7 - ROH")
cli::cli_alert_info("Setup environment")
suppressMessages({
  library(tidyverse)
  library(here)
  library(patchwork)

  source(here("R", "setup.R"))
  source(here("R", "plink.R"))
  source(here("R", "sample.R"))
})

## Functions -----
length_breaks <- c(0, 500, 1e3, 2e3, 4e3, 8e3, 16e3, Inf)
length_labels <- c("0-.5", ".5-1", "1-2", "2-4", "4-8", "8-16", ">16")

subgroup_subset <- c(
  "Africa", "West Eurasia", "East Asia", "North America",
  "Amazon", "West Andes", "Southeast America", "Southern America"
)

summarise_roh <- function(dat, ...) {
  dat %>%
    group_by(IID, ...) %>%
    summarise(total_length = sum(KB)/1e3, rohs = n(), .groups =  "drop") %>%
    ungroup()
}

plot_roh <- function(dat, shape) {
  studies <- c("Present study", "SGDP", "Crawford2017", "Fuente2018", "Raghavan2015")
  ## Define default setup
  plot_pointrange <- list(
    stat_summary(
      fun.data = median_hilow, geom = "linerange",
      position = position_dodge(width = .8)
    ),
    stat_summary(
      fun = median, geom = "point", size = 1.5,
      position = position_dodge(width = .8)
    )
  )
  plot_scatter <- geom_point(size = 1.5)

  plot_style <- list(
    scale_shape_manual("Study", values = set_names(1:5 + 20, studies)),
    scale_color_manual("Region", values = subgroups_pal2),
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0),
      strip.background = element_blank()
    )
  )

  roh_pop <- summarise_roh(dat, {{shape}}, pop, subgroup) %>%
    mutate(study = parse_factor(study, studies)) %>%
    filter(subgroup == "Amazon") %>%
    ggplot(aes(total_length, pop, shape = {{shape}}, color = subgroup)) +
    theme_classic(8, "Helvetica") +
    plot_pointrange +
    plot_style +
    guides(shape = "none", color = "none") +
    # facet_grid(subgroup~., scales = "free", space = "free") +
    labs(y = "Population", x = "Total ROH length (Mbp)")

  roh_length <- summarise_roh(dat, {{shape}}, subgroup, length) %>%
    mutate(study = parse_factor(study, studies)) %>%
    ggplot(aes(length, total_length, color = subgroup, shape = {{shape}})) +
    plot_pointrange +
    plot_style +
    guides(shape = "none", color = "none") +
    labs(x = "ROH length category (Mbp)", y = "Total ROH length (Mbp)")

  roh_scatter <- summarise_roh(dat, {{shape}}, subgroup) %>%
    mutate(study = parse_factor(study, studies)) %>%
    ggplot(aes(total_length, rohs, color = subgroup, shape = {{shape}})) +
    plot_scatter +
    plot_style +
    guides(
      shape = ggplot2::guide_legend(title.position = "top", nrow = 2),
      color = ggplot2::guide_legend(title.position = "top", nrow = 3)
    ) +
    labs(x = "Total ROH length (Mbp)", y = "Total number of ROHs")

  (roh_length + roh_scatter + guide_area() + roh_pop) +
    plot_layout(guides = "collect", design = "14\n23") +
    plot_annotation(tag_levels = "A") &
    theme(legend.direction = "horizontal")
}

## Compute ROH -----
if (!file.exists(here("out", "roh", "DataB.hom"))) {
  cli::cli_alert_info("Computing ROH data")

  dir.create(here("out", "roh"), recursive = TRUE, showWarnings = FALSE)
  res <- exec_plink(c(
    "-bfile", here("out/plink/DataB"),
    "--homozyg",
    "--homozyg-kb", "500",
    "--homozyg-gap", "100",
    "--homozyg-density", "50",
    "--homozyg-window-snp", "50",
    "--homozyg-window-het", "1",
    "--homozyg-window-missing", "25",
    "--homozyg-window-threshold", ".05",
    "--out", "out/roh/DataB"
  ))
  if (res != 0) {
    stop("Unable to run plink to generate ROH data")
  }
  cli::cli_alert_success("Completed ROH data calculation")
}

## Load data -----
cli::cli_alert_info("Loading ROH data")

roh <- here("out", "roh", "DataB.hom") %>%
  read_table2(col_types = "ccdicciididdd") %>%
  left_join(read_sample_info(), by = c(IID = "sample")) %>%
  filter(subgroup %in% subgroup_subset) %>%
  mutate(
    length = cut(KB, length_breaks, length_labels)
  )

## Plots -------
# ## A) ROH by age
# cli::cli_alert_info("Generating ROH plots by sample age")
# save_plot(
#   here("figs", "sup-fig", "sup-fig_roh-age"),
#   plot_roh(roh, age), width = 8, height = 8
# )

## B) ROH Study
cli::cli_alert_info("Generating ROH plots by sample study")
save_plot(
  here("figs", "sup-fig", "sup-fig3_roh"),
  plot_roh(filter(roh, age == "Contemporan"), study),
  width = 8, height = 6
)
cli::cli_alert_success("Done!")
