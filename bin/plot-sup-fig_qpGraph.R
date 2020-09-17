cli::cli_h1("Plotting Sup. Figure 8-10 - qpGraph")
cli::cli_alert_info("Setup environment")

suppressMessages({
  library(tidyverse)
  library(here)
  library(tidygraph)
  library(ggraph)
  library(patchwork)
  library(ggbeeswarm)
  library(admixtools)

  source(here("R", "setup.R"))
  source(here("R", "plink.R"))
  source(here("R", "sample.R"))
})

dir.create(here("figs", "sup-fig"), showWarnings = FALSE)

## Functions -----
admixtree_layout <- function(edge) {
  edge <- mutate(edge, weight = if_else(type != "admix", 100 * weight, weight))
  node <- mutate(node, region = name)
  extd <- tibble(
    name = setdiff(edge$to, node$name),
    region = edge$from[match(name, edge$to)]
  )
  layout <-
    tbl_graph(nodes = bind_rows(extd, node), edges = edge) %>%
    create_layout("tree")
  admix_node <- filter(edge, type == "admix") %>%
    mutate(fy = layout$y[match(from, layout$name)]) %>%
    group_by(to) %>%
    summarise(y = min(fy) - 1)
  nodes <- admix_node$to
  new_y <- admix_node$y
  while(length(nodes) > 0) {
    layout$y[match(nodes, layout$name)] <- new_y
    nodes <- edge$to[match(nodes, edge$from)]
    new_y <- new_y - 1
    nodes <- nodes[!is.na(nodes)]
    new_y <- new_y[!is.na(nodes)]
  }
  return(layout)
}

plot_admixtree <- function(layout, ...) {
  ggraph(layout, color = "darkgray") + 
    geom_edge_link(
      aes(linetype = type, ...),
      angle_calc = 'along', label_dodge = unit(2.5, 'mm'), label_size = 2
    ) +
    geom_node_label(aes(
      label = name, color = region,
      filter = !str_starts(name, "SA") & name != "Root"
    ), size = 2) +
    scale_edge_linetype_manual(
      "Edge",
      values = c(edge = "solid", admix = "dashed"),
      labels = c(admix = "Admixture", edge = "Normal"),
      limits = c("edge", "admix")
    ) +
    scale_color_manual(
      values = subgroups_pal, guide = "none", na.value = "darkgray"
    ) +
    theme(legend.position = c(0, 1), legend.justification = c(0, 1)) +
    coord_cartesian(clip = "off")
}

add_leaves <- function(model, leaves) {
  tibble(from = names(leaves), to = as.character(leaves)) %>%
    bind_rows(model) %>%
    as_tbl_graph()
}

## Setup -----
cli::cli_alert_info("Load data")

data_prefix <- here("out", "plink", "DataB")
data_f2_dir <- here("out", "qpF2")

fam <- read_fam(paste0(data_prefix, ".fam")) %>%
  left_join(read_sample_info()) %>%
  filter(
    age == "Contemporan",
    group %in% c("NAT", "OCE") | pop %in% c("Mbuti", "Han")
  )

## Define Models --------
node <- tibble(
  name = c(
    "Root", "SAB", "SA1", "SA2", "SA3", "Out",
    "WestAndes", "Amazon", "South", "SouthEast"
  ),
  region = name
)

models <- list(
  model1 = tribble(
    ~ from, ~ to, ~ type,
    "Root", "Out", "edge",
    "Root", "SAB", "edge",
    "SAB", "SA1", "edge",
    "SAB", "SA2", "edge",
    "SA1", "South", "edge",
    "SA1", "WestAndes", "edge",
    "SA2", "Amazon", "edge",
    "SA2", "SouthEast", "edge",
  ),
  model2 = tribble(
    ~ from, ~ to, ~type,
    "Root", "Out", "edge",
    "Root", "SA1", "edge",
    "SA1", "South", "edge",
    "SA1", "SA2", "edge",
    "SA2", "SA3", "edge",
    "SA2", "Amazon", "edge",
    "SA3", "WestAndes", "edge",
    "SA3", "SouthEast", "edge",
  ),
  model3 = tribble(
    ~ from, ~ to, ~type,
    "Root", "Out", "edge",
    "Root", "SA1", "edge",
    "SA1", "SA2", "edge",
    "SA1", "WestAndes", "edge",
    "SA2", "SA3", "edge",
    "SA2", "South", "edge",
    "SA3", "Amazon", "edge",
    "SA3", "SouthEast", "edge",
  )
)

## Pre-Compute model likelihood ----
model_dat_file <- here("out", "model_qpgraph.rds")
if (!file.exists(model_dat_file)) {
  cli::cli_alert_info("Running qpGraph")
  suppressMessages({
    library(furrr)
    plan(multisession, workers = 4)
  })
  ## Define populations to investigate
  pop_subset <- unique(fam$pop)
  pop_group <- distinct(fam, pop, subgroup) %>% with(split(pop, subgroup))
  pop_combn <- list(
    Out = c("Mbuti", "Han", "Pima", "Mayan"),
    WestAndes = pop_group[["West Andes"]],
    Amazon = pop_group[["Amazon"]],
    South = pop_group[["South"]],
    SouthEast = pop_group[["South East"]]
  ) %>%
  cross()

  ## Pre-compute F2
  if (!dir.exists(data_f2_dir)) {
    cli::cli_alert_info("Running F2 computation")
    extract_f2(data_prefix, data_f2_dir, pops = pop_subset, overwrite = TRUE)
  }
  data_f2 <- f2_from_precomp(data_f2_dir)

  cli::cli_alert_info("Fitting qpGraph models")
  dat <- list(combn = pop_combn, model_name = names(models)) %>%
    cross_df() %>%
    mutate(
      model = map(model_name, ~ models[[.x]]) %>% map2(combn, add_leaves),
      combn = map(combn, as_tibble),
      graph = future_map(
        model, qpgraph, f2_blocks = data_f2, return_f4 = TRUE, .progress = TRUE
      )
    ) %>%
    unnest(combn)
  saveRDS(dat, model_dat_file)
}
dat <- readRDS(model_dat_file)

## Plots ----
cli::cli_alert_info("Generating Sup Figure 8 - Models overview")
model_plots <- models %>%
  map(tbl_graph, nodes = node) %>%
  map(create_layout, layout = "tree") %>%
  map(plot_admixtree) %>%
  wrap_plots() +
  plot_annotation(tag_prefix = "Model ", tag_levels = "1") +
  plot_layout(guides = "collect")
save_plot(
  here("figs", "sup-fig", "sup-fig8_qpgraph_models"),
  model_plots, width = 8, height = 3
)

cli::cli_alert_info("Generating Sup Figure 9 - Top 10 Models")
top_models <- dat %>%
  mutate(
    score = map_dbl(graph, "worst_residual"),
    edges = map(graph, "edges"),
    name = sprintf("%s - Worst F4: %.2f", model_name, score)
  ) %>%
  slice_min(score, n = 10)

top_models_plot <- top_models$edges %>%
  map(admixtree_layout) %>%
  map(plot_admixtree, label = round(weight, 2)) %>%
  map2(top_models$name, ~ .x + ggtitle("", subtitle = .y)) %>%
  wrap_plots(ncol = 5) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "1", tag_prefix = "#") &
  theme(legend.position = "bottom")

save_plot(
  here("figs", "sup-fig", "sup-fig9_qpgraph_top_models"),
  top_models_plot, width = 10, height = 8
)

cli::cli_alert_info("Generating Sup Figure 10 - Models score")
plot_score <- function(dat, score = "score") {
  dat %>%
  mutate(Out = parse_factor(Out, c("Mbuti", "Han", "Mayan", "Pima"))) %>%
  mutate(score = map_dbl(graph, score)) %>%
  ggplot(aes(model_name, abs(score))) +
  geom_quasirandom(aes(color = model_name), show.legend = FALSE) +
  stat_summary(fun.data = median_hilow) +
  facet_grid(Out ~ ., labeller = label_both) +
  scale_x_discrete("") +
  scale_color_brewer(palette = "Set1") +
  coord_flip() +
  theme(strip.background.y = element_blank(), strip.text.y = element_blank())
}

model_lik <- plot_score(dat) +
  scale_y_log10("Likelihood score")
model_f4  <- plot_score(dat, "worst_residual") +
  scale_y_log10("Worst Absolute F4") +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )

save_plot(
  here("figs", "sup-fig", "sup-fig10_qpgraph_model_likelihood"),
  model_lik + model_f4, width = 8, height = 4
)

cli::cli_alert_success("Done!")