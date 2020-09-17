set.seed(551991)

options(dplyr.summarise.inform = FALSE)

## Wrapper to save both PDF and TIFF
save_plot <- function(prefix, plot, ...) {
  dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
  
  ggsave(paste0(prefix, ".pdf"), plot, ..., useDingbats = FALSE)
  ggsave(
    paste0(prefix, ".tiff"), plot, ...,
    dpi = 300, compression = "lzw", type = "cairo"
  )
  ggsave(paste0(prefix, ".png"), plot, ..., dpi = 150)
}

## Setup environment
ggplot2::theme_set(ggplot2::theme_classic(10, "Helvetica"))

cbpal <- c(
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f",
  "#ffffff", "#333333"
)

groups <- c("AFR", "EUR", "EAS", "NAT", "OCE", "CAS", "SAS")

subgroups <- c(
  "North America", "Central America", "Amazon",
  "South East", "West Andes", "South", "Eskimo"
)
subgroups_pal <-
  c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69")
names(subgroups_pal) <- subgroups

subgroups_pal2 <-set_names(
  cbpal[1:10],
  c(subgroups, "Africa", "West Eurasia", "East Asia")
)

subgroup_full <- c(
  "Africa", "West Eurasia", "South West Asia", "East Asia", "Eskimo",
  "North America", "Central America", "Amazon", "West Andes",
  "South East", "South", "Central Asia & Siberia", "Oceania",
  "Papua New Guinea"
)

subgroups2 <- c(
  "Africa", "West Eurasia", "South West Asia", "East Asia",
  subgroups, paste0("Ancient\n", subgroups),
  "Central Asia & Siberia", "Oceania", "Papua New Guinea"
  ) %>%
  str_replace("\n(South|Eskimo)", " \\1")
