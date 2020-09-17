cli::cli_h1("Plotting Sup. Figure 2 - Treemix")
## Setup environment -----
cli::cli_alert_info("Setup environment")
suppressMessages({
  library(tidyverse)
  library(here)
  library(ggtree)
  library(ape)
  library(patchwork)

  source(here("R", "setup.R"))
  source(here("R", "sample.R"))
  source(here("R", "plink.R"))
  source(here("R", "dstat.R"))
})

dir.create(here("figs", "sup-fig"), recursive = TRUE, showWarnings = FALSE)

## Function ------
exec_treemix <- function(
  input, output, block = 50, root = "Mbuti", migration = NULL,
  seed = 13121011, force = FALSE, stdout = stdout()
) {
  suppressMessages(require(sys))
  treemix <- c(
    "docker", "run", "--rm", "-w", here(), "-v", paste0(here(), ":", here()),
    "quay.io/biocontainers/treemix:1.13--h4bb999f_2", "treemix"
  )
  if (!force && file.exists(paste0(output, ".treeout.gz")))
    return(0L)

  arg <- c("-i", input, "-o", output, "-k", block, "-root", root)
  if (!is.null(seed))
    arg <- c(arg, "-seed", seed)
  if (!is.null(migration))
    arg <- c(arg, "-m", migration)
  sys::exec_wait(treemix[1], c(treemix[-1], arg), std_out = stdout)
}

add_conf <- function(tree, boots) {
  treedata <- fortify(tree)
  treedata$conf <- NA_real_
  treedata$conf[is.na(treedata$label)] <- ape::prop.clades(tree, boots)
  treedata
}

tree_crop <- function(tree, crop = .05) {
  tree$edge.length <- pmin(tree$edge.length, crop)
  tree
}

read_migration <- function(tree, lines) {
  node_df <- fortify(tree) %>% select(node, x, y)

  extract_node <- function(node) {
    if (str_detect(node, "\\(")) {
      nodes <- str_extract_all(node, "[A-Za-z]+")[[1]]
      nodes <- nodes[nodes %in% tree$tip.label]
      return(MRCA(tree, nodes))
    }
    return(match(str_extract(node, "[A-Za-z]+"), tree$tip.label))
  }
  readr::read_delim(
    lines[-1], delim = " ",
    col_names = c("rate", "A", "B", "C", "source", "dest")
    ) %>%
    mutate(
      source_node = map_int(source, extract_node),
      dest_node = map_int(dest, extract_node)
    ) %>%
    left_join(node_df, by = c(source_node = "node")) %>%
    left_join(
      rename(node_df, xend = "x", yend = "y"),
      by = c(dest_node = "node")
    )
}

## Load data -----
cli::cli_alert_info("Loading data")
dat <- read_plink(here("out", "plink", "DataB"))
fam <- left_join(dat$fam, read_sample_info(), by = c("pop", "sample"))

## Write treemix -----
file <- here("out", "treemix", "DataB_fq01.treemix.gz")

if (!file.exists(file)) {
  cli::cli_alert_warning("Unable to find treemix dataset")
  cli::cli_alert_info("Running treemix")

  fam_subset <- fam %>%
    filter(age == "Contemporan", group == "NAT" | pop == "Mbuti") %>%
    mutate(idx = match(sample, rownames(dat$bed)))
  samples <- fam_subset$idx
  populations <- with(fam_subset, split(idx, pop))
  common_vars <- which(allele_freq(samples, dat$bed) > 0.01)

  cli::cli_alert_info(c(
    "Processing ",
    "{length(common_vars)} variant{?s} among ",
    "{length(samples)} sample{?s} within ",
    "{length(populations)} population{?s}."
  ))

  ac <- map(populations, allele_count, mtx = dat$bed[, common_vars])
  an <- map(populations, allele_total, mtx = dat$bed[, common_vars])
  treemix <- names(populations) %>%
    set_names(.) %>%
    map(~ sprintf("%d,%d", ac[[.x]], an[[.x]] - ac[[.x]])) %>%
    bind_cols()

  ## Save
  dir.create(here("out", "treemix"), showWarnings = FALSE, recursive = TRUE)
  write_tsv(treemix, file)
  cli::cli_alert_success("Saved treemix data at `{file}`")
}

## Run Treemix ------
best <- here("out", "treemix", "DataB_fq01")
boot <- here("out", "treemix", "DataB_fq01_boot")

if (!file.exists(paste0(best, ".treeout.gz"))) {
  cli::cli_alert_warning("Generating treemix ML tree")
  exec_treemix(input = file, output = best, migration = 5, block = 50)
  cli::cli_alert_success("Complete treemix run")
}

if (!file.exists(paste0(boot, ".treeout.gz"))) {
  cli::cli_alert_warning("Generating treemix ML tree bootstraps")
  suppressMessages({
    library(furrr)
    plan(multisession, workers = 4)
  })

  files <- paste0(boot, "-", 1:500)
  reslt <- future_map_int(
    files, exec_treemix, input = file, block = 50, seed = NULL, stdout = FALSE, force = TRUE,
    .progress = TRUE
  )
  if (any(reslt != 0)) {
    failed <- sum(reslt != 0)
    cli::cli_alert_danger("{sum(reslt != 0)} bootstrap runs failed!")
    stop(str_glue("A total of {failed} bootstrap run failed!"))
  }

  map(files, paste0, ".treeout.gz") %>%
    map(read_lines) %>%
    reduce(c) %>%
    write_lines(paste0(boot, ".treeout.gz"))
  cli::cli_alert_success("Complete treemix bootstrap")
}

## read sample info ------
cli::cli_alert_info("Generating treemix plots")
pop_subgroup <- fam %>% distinct(pop, subgroup) %>% with(split(pop, subgroup))

treemix_lines <- read_lines(paste0(best, ".treeout.gz"))
treemix_tree <- read.tree(text = treemix_lines[1])
treemix_migr <- read_migration(treemix_tree, treemix_lines)
treemix_boot <- read.tree(paste0(boot, ".treeout.gz"))

treemix_best <- treemix_tree %>%
  groupOTU(pop_subgroup) %>%
  add_conf(treemix_boot) %>%
  mutate(group = if_else(group == "0", NA_character_, as.character(group))) %>%
  ggtree(size = NA) +
  geom_tree(aes(color = group)) +
  geom_tiplab(size = 2.5) +
  geom_nodelab(
    size = 2, hjust = 1.5, vjust = 0,
    aes(label = if_else(conf / 5 >= 75, conf / 5, NA_real_))
  ) +
  geom_curve(
    data = treemix_migr,
    arrow = arrow(length = unit(0.02, "npc")),
    aes(xend = xend, yend = yend, alpha = rate)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, .15))) +
  scale_color_manual(
    "Region", values = subgroups_pal2, limits = c("Africa", subgroups[-7]),
    na.value = "darkgray"
  ) +
  scale_alpha("Migration rate") +
  theme_tree2() +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key.height = unit(.8, "lines")
  )
treemix_consensus <- consensus(treemix_boot, p = 0.75) %>%
  groupOTU(pop_subgroup) %>%
  add_conf(treemix_boot) %>%
  mutate(group = if_else(group == "0", NA_character_, as.character(group))) %>%
  ggtree(size = NA) +
  geom_tree(aes(color = group)) +
  geom_tiplab(size = 2.5) +
  geom_nodelab(
    size = 2, hjust = 1.5, vjust = 0,
    aes(label = if_else(conf / 5 > 75, conf / 5, NA_real_))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_color_manual(
    "Region",
    values = subgroups_pal2,
    limits = c("Africa", subgroups[-length(subgroups)]),
    na.value = "darkgray",
  ) +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key.height = unit(.8, "lines")
  )

## Compile as a single figure
save_plot(
  here("figs", "sup-fig", "sup-fig2_treemix"),
  treemix_best / treemix_consensus,
  width = 6, height = 8
)

cli::cli_alert_success("Done!")