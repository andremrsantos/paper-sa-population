library(tidyverse)
library(here)
library(patchwork)

source(here("R", "setup.R"))
source(here("R", "sample.R"))

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
    dplyr::summarise(total_length = sum(KB)/1e3, rohs = n()) 
}

legend_style <- function(at = c(1, 1), color = 1, shape = 2) {
  list(
    ggplot2::guides(
      color = guide_legend(order = color),
      shape = guide_legend(order = shape)
    ),
    ggplot2::theme(
      legend.position = at,
      legend.justification = at,
      legend.background = element_blank(),
      legend.box = "horizontal",
      legend.title = element_text(size = 5), 
      legend.text = element_text(size = 5),
      legnd.key.size = unit(.75, "lines")
    )
  )
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

#   "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
#   "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"
# )

roh <- here("out", "roh", "DataB.hom") %>%
  read_table2(col_types = "ccdicciididdd") %>%
  left_join(read_sample_info(), by = c(IID = "sample")) %>%
  filter(subgroup %in% subgroup_subset) %>%
  mutate(
    length = cut(KB, length_breaks, length_labels)
  )

## Plots -------
## A) ROH by age
roh_age <- summarise_roh(roh, age, subgroup, length) %>%
  ggplot(aes(length,total_length, color = subgroup, shape = age)) +
  stat_summary(fun.data = median_hilow, position = position_dodge(width = .8)) +
  scale_shape_manual(values = c(21, 22)) +
  scale_color_manual(values = subgroup_palette) +
  labs(
    x = "ROH length category (Mbp)",
    y = "Total ROH length (Mbp)",
    shape = "Sample Age",
    color = "Region"
  ) +
  legend_style(color = 2, shape = 1)

roh_age_scatter <- summarise_roh(roh, age, subgroup) %>%
  ggplot2::ggplot(aes(total_length, rohs, color = subgroup, shape = age)) +
  ggplot2::geom_point() +
  ggplot2::scale_shape_manual(values = c(21, 22)) +
  ggplot2::scale_color_manual(values = subgroup_palette) +
  ggplot2::labs(
    x = "Total ROH length (Mbp)",
    y = "Total number of ROHs",
    shape = "Sample Age",
    color = "Region"
  ) +
  legend_style(c(0, 1))

ggsave(
  here("figs", "sup-fig", "sup-fig_roh-age.pdf"),
  (roh_age + roh_age_scatter) + plot_annotation(tag_levels = "A"),
  width = 8, height = 4, useDingbats = FALSE
)

## B) ROH Study
roh_study <- roh %>%
  filter(age == "Contemporan") %>%
  summarise_roh(study, subgroup, length) %>%
  ggplot(aes(length, total_length, color = subgroup, shape = study)) +
  stat_summary(fun.data = median_hilow, position = position_dodge(width = .8)) +
  scale_shape_manual(values = 1:5 + 20) +
  scale_color_manual(values = subgroup_palette) +
  labs(
    x = "ROH length category (Mbp)",
    y = "Total ROH length (Mbp)",
    shape = "Study",
    color = "Region"
  ) +
  legend_style(color = 2, shape = 1)

roh_study_scatter <- roh %>%
  filter(age == "Contemporan") %>%
  summarise_roh(study, subgroup) %>%
  ggplot(aes(total_length, rohs, color = subgroup, shape = study)) +
  geom_point() +
  scale_shape_manual(values = 1:5 + 20) +
  scale_color_manual(values = subgroup_palette) +
  labs(
    x = "Total ROH length (Mbp)",
    y = "Total number of ROHs",
    shape = "Sample Age",
    color = "Region"
  ) +
  legend_style(c(0, 1))

ggsave(
  here("figs", "sup-fig", "sup-fig_roh-study.pdf"),
  (roh_study + roh_study_scatter) + plot_annotation(tag_levels = "A"),
  width = 8, height = 4, useDingbats = FALSE
)
