suppressMessages({
  library(tidyverse)
  library(here)
  library(patchwork)

  source(here("R", "setup.R"))
  source(here("R", "plink.R"))
  source(here("R", "sample.R"))
})

dir.create(here("figs", "sup-fig"), recursive = TRUE, showWarnings = FALSE)

## Compute ROH -----
if (!file.exists(here("out", "roh", "DataB.hom"))) {
  dir.create(here("out", "roh"), recursive = TRUE, showWarnings = FALSE)
  exec_plink(c(
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
}

## Functions -----
summarise_roh <- function(dat, ...) {
  dat %>%
    dplyr::group_by(IID, ...) %>%
    dplyr::summarise(total_length = sum(KB)/1e3, rohs = n(), .groups =  "drop") 
}

# legend_style <- function(at = c(1, 1), color = 1, shape = 2) {
#   list(
#     ggplot2::guides(
#       color = guide_legend(order = color),
#       shape = guide_legend(order = shape)
#     ),
#     ggplot2::theme(
#       legend.position = at,
#       legend.justification = at,
#       legend.background = element_blank(),
#       legend.box = "horizontal",
#       legend.title = element_text(size = 5), 
#       legend.text = element_text(size = 5),
#       legend.key.height = unit(.75, "lines")
#     )
#   )
# }

plot_roh <- function(dat, shape) {
  ## Define default setup
  plot_pointrange <- list(
    stat_summary(fun.data = median_hilow, position = position_dodge(width = .8), geom = "linerange"),
    stat_summary(fun = median, position = position_dodge(width = .8), geom = "point", size = 1.5)
  )
  plot_scatter <- geom_point(size = 1.5)

  plot_style <- list(
    scale_shape_manual("Sample Age", values = 1:5 + 20),
    scale_color_manual("Region", values = subgroup_palette),
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0),
      strip.background = element_blank()
    )
  )

  roh_pop <- summarise_roh(roh, {{shape}}, pop, subgroup) %>%
    ggplot(aes(total_length, pop, shape = {{shape}}, color = subgroup)) +
    theme_classic(8, "Helvetica") +
    plot_pointrange +
    plot_style +
    facet_grid(subgroup~., scales = "free", space = "free") +
    labs(y = "Population", x = "Total ROH length (Mbp)")
  
  roh_length <- summarise_roh(roh, {{shape}}, subgroup, length) %>%
    ggplot(aes(length, total_length, color = subgroup, shape = {{shape}})) +
    plot_pointrange +
    plot_style +
    guides(shape = "none", color = "none") +
    labs(x = "ROH length category (Mbp)", y = "Total ROH length (Mbp)")

  roh_scatter <- summarise_roh(roh, {{shape}}, subgroup) %>%
    ggplot(aes(total_length, rohs, color = subgroup, shape = {{shape}})) +
    plot_scatter +
    plot_style +
    guides(shape = "none", color = "none") +
    labs(x = "Total ROH length (Mbp)", y = "Total number of ROHs")
  
  ((roh_length / roh_scatter) | roh_pop) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A")
}

## Load data -----
length_breaks <- c(0, 500, 1e3, 2e3, 4e3, 8e3, 16e3, Inf)
length_labels <- c("0-.5", ".5-1", "1-2", "2-4", "4-8", "8-16", ">16")

subgroup_subset <- c(
  "Africa", "West Eurasia", "East Asia", "North America",
  "Amazon", "West Andes", "South East", "South"
)
subgroup_palette <- c(
  subgroups_pal,
  set_names(setdiff(cbpal, subgroups_pal)[1:3], subgroup_subset[1:3])
)

roh <- here("out", "roh", "DataB.hom") %>%
  read_table2(col_types = "ccdicciididdd") %>%
  left_join(read_sample_info(), by = c(IID = "sample")) %>%
  filter(subgroup %in% subgroup_subset) %>%
  mutate(
    length = cut(KB, length_breaks, length_labels)
  )

## Plots -------
## A) ROH by age
save_plot(
  here("figs", "sup-fig", "sup-fig_roh-age"),
  plot_roh(roh, age), width = 8, height = 10
)

## B) ROH Study
save_plot(
  here("figs", "sup-fig", "sup-fig_roh-study"),
  plot_roh(filter(roh, age == "Contemporan"), age),
  width = 8, height = 10
)
