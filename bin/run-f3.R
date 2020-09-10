library(tidyverse)
library(here)

source(here("R", "setup.R"))
source(here("R", "sample.R"))
source(here("R", "plink.R"))
source(here("R", "dstat.R"))
source(here("R", "map.R"))

dir.create(here("figs", "sup-fig"), recursive = TRUE, showWarnings = FALSE)

## Load data -----
plink <- read_plink(here("out", "plink", "DataB"))

dat <- plink$fam %>%
  left_join(read_sample_info()) %>%
  filter(group %in% c("NAT", "AFR"), subgroup != "Africa" | pop == "Mbuti")

pop_info <- dat %>%
  group_by(pop, subgroup, age) %>%
  summarise(across(c(long, lat), mean))

pop_sample <- dat %>%
  with(c(split(sample, pop), split(sample, subgroup))) %>%
  discard(~ length(.x) == 0)
pop_count <- compute_count(plink, pop_sample)
pop_total <- compute_total(plink, pop_sample)

## Compute F3 -----
south  <- unique(filter(dat, subgroup == "South", age == "Contemporan")$pop)
amazon <- unique(filter(dat, subgroup == "Amazon", age == "Contemporan")$pop)
wandes <- unique(filter(dat, subgroup == "West Andes", age == "Contemporan")$pop)
native <- unique(filter(dat, group == "NAT", age == "Contemporan")$pop)

# f3_admix <- list(
#   south_amazon = f3_stat(native, south, amazon, pop_count, pop_total),
#   south_wandes = f3_stat(native, south, wandes, pop_count, pop_total),
#   amazon_wandes = f3_stat(native, amazon, wandes, pop_count, pop_total)
# )

native <- unique(filter(dat, group == "NAT")$pop)

f3_outgroup_file <- here("out", "f3-outgroup.rds")
if (file.exists(f3_outgroup_file)) {
  f3_outgroup <- readRDS(f3_outgroup_file)
} else {
  f3_outgroup <- f3_stat("Mbuti", native, native, pop_count, pop_total)
  saveRDS(f3_outgroup, f3_outgroup_file)
}

## Plotting -----
f3_label_pop <- pop_info %>% 
  filter(age == "Ancient", pop %in% f3_outgroup$y) %>%
  mutate(y = paste0("Y=", pop))

f3_map <- f3_outgroup %>%
  filter(y %in% f3_label_pop$pop) %>%
  left_join(pop_info, by = c(x = "pop")) %>%
  group_by(x, y, subgroup, age, lat, long) %>%
  summarise(f3 = mean(zv)) %>%
  group_by(y) %>%
  mutate(f3 = rank(f3), y = paste0("Y=", y)) %>%
  ggplot() +
  base_america_map() +
  geom_point(aes(long, lat, color = f3, shape = age)) +
  geom_point(aes(long, lat), shape = 3, color = "firebrick", data = f3_label_pop) +
  scale_color_viridis_c() +
  facet_wrap(~ y) +
  labs(x = "", y = "", title = "F3(Mbuti, X, Y)") +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
save_plot(
  here::here("figs", "sup-fig", "sup-fig_f3-map"),
  f3_map, height = 10, width = 8
)
