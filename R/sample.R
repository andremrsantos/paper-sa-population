library(tidyverse)
library(here)

source(here("R", "setup.R"))

read_sample_info <- function() {
  read_tsv(here("data", "sample-info.tsv"), col_types = cols()) %>%
  filter(study != "Skoglund2015") %>%
  mutate(
    age = parse_factor(sample_age, c("Contemporan", "Ancient")),
    group = parse_factor(group, groups),
    subgroup = parse_factor(subgroup, subgroup_full),
    subgroup2 = if_else(
      age == "Ancient",
      paste0("Ancient\n", subgroup),
      as.character(subgroup)
    ) %>%
      str_replace("\n(Eskimo)", " \\1") %>%
      parse_factor(subgroups2)
  ) %>%
  select(-sample_age)
}
