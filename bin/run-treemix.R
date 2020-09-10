library(tidyverse)
library(here)
library(matrixStats)
## Read treeout
# BiocManager::install("ggtree")
library(ggtree)
library(ape)
library(patchwork)

source(here("R", "setup.R"))
source(here("R", "sample.R"))
source(here("R", "plink.R"))
source(here("R", "dstat.R"))

dir.create(here("figs", "sup-fig"), recursive = TRUE, showWarnings = FALSE)

## Function ------
exec_treemix <- function(
  input, output, block = 50, root = "Mbuti", migration = NULL
  ) {
  require(sys)
  treemix <- c(
    "docker", "run", "--rm", "-w", here(), "-v", paste0(here(), ":", here()),
    "quay.io/biocontainers/treemix:1.13--h4bb999f_2", "treemix"
  )
  
  arg <- c("-i", input, "-o", output, "-k", block, "-root", root)
  if (!is.null(migration))
    arg <- c(arg, "-m", migration)
  sys::exec_wait(treemix[1], c(treemix[-1], arg), std_out = FALSE)
}

add_conf <- function(tree, boots) {
  treedata <- fortify(tree)
  treedata$conf <- NA_real_
  treedata$conf[is.na(treedata$label)] <- ape::prop.clades(tree, boots)
  treedata
}

read_migration <- function(lines) {
  tree <- read.tree(text = lines[1])
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
    left_join(rename(node_df, dest_node = "node", xend = "x", yend = "y"))
}

## Load data -----
dat <- read_plink(here("out", "plink", "DataB"))
fam <- left_join(dat$fam, read_sample_info())

## Write treemix -----
file <- here("out", "treemix", "DataB_fq01.treemix.gz")

if (!file.exists(file)) {
  idv <- filter(fam, age == "Contemporan", group == "NAT" | pop == "Mbuti") %>%
    pull(sample) %>%
    match(rownames(dat$bed))
  print("No Samples: ")
  print(length(idv))
  
  pop <- filter(fam, age == "Contemporan", group == "NAT" | pop == "Mbuti") %>%
    with(split(sample, pop)) %>%
    map(match, rownames(dat$bed))
  print("No Populations: ")
  print(length(pop))
  
  ac <- allele_count
  vr <- which(colMeans2(dat$bed, rows = idv, na.rm = TRUE)/2 > 0.01)
  an <- purrr::map(pop, ~ 2 * length(.x) - colCounts(dat$bed, cols = vr, rows = .x, value = NA_integer_))
  a1 <- purrr::map(pop, ~ colSums2(dat$bed, rows = .x, na.rm = TRUE))
  treemix <- names(pop) %>%
    purrr::set_names(.) %>%
    purrr::map(~ paste0(a1[[.x]], ",", an[[.x]] - a1[[.x]])) %>%
    dplyr::bind_cols()
  ## Save
  dir.create(here("out", "treemix"), showWarnings = FALSE, recursive = TRUE)
  write_tsv(treemix, file)
}

## Run Treemix
best <- here("out", "treemix", "DataB_fq01")
boot <- here("out", "treemix", "DataB_fq01_boot")

if (!file.exists(paste0(best, ".treeout.gz"))) {
  exec_treemix(input = file, output = best, migration = 5)
}

if (!file.exists(paste0(boot, ".treeout.gz"))) {
  library(furrr)
  plan(multisession, workers = 4)
  
  files <- paste0(boot, '-', 1:5)
  reslt <- future_map_int(files, exec_treemix, input = file, .progress = TRUE)
  if (any(reslt != 0)) {
    failed <- sum(reslt != 0)
    stop(str_glue("A total of {failed} bootstrap run failed!"))
  }
  map(files, paste0, ".treeout.gz") %>%
    map(read_lines) %>%
    reduce(c) %>%
    write_lines(paste0(boot, ".treeout.gz"))
}

## read sample info ------
pop_subgroup <- fam %>% distinct(pop, subgroup) %>% with(split(pop, subgroup))

treemix_lines <- read_lines(paste0(best, ".treeout.gz"))
treemix_tree <- read.tree(text = treemix_lines[1])
treemix_boot <- read.tree(paste0(boot, ".treeout.gz"))
treemix_migr <- read_migration(treemix_lines)

treemix_best <-
  treemix_tree %>%
  groupOTU(pop_subgroup) %>%
  add_conf(treemix_boot) %>%
  ggtree(size = NA) +
  geom_tree(layout = "slanted", aes(color = group)) +
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
  scale_color_manual("Region", values = subgroups_pal2) +
  scale_alpha("Migration rate") +
  theme_tree2() +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key.height = unit(.8, "lines")
  )
# ggsave(
#   here("figs", "sup-fig_treemix.pdf"),
#   treemix_best, width = 4.5, height = 4, useDingbats = FALSE
# )

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
# ggsave(
#   here("figs", "sup-fig_treemix-consensus.pdf"),
#   treemix_consensus, width = 4.5, height = 4, useDingbats = FALSE
# )

## Compile as a single figure
save_plot(
  here("figs", "sup-fig", "sup-fig_treemix-merged"),
  treemix_best / treemix_consensus,
  width = 6, height = 8
)



