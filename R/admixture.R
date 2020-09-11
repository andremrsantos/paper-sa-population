library(tidyverse)

read_q <- function(file, fam, group) {
  dat <- file %>%
    read_table(col_names = FALSE, col_types = cols()) %>%
    mutate(group = pull(fam, {{group}}))
  
  k_range <- seq_len(ncol(dat) - 1)
  k_means <- dat %>% group_by(group) %>% summarise(across(.fns = mean))
  k_order <- apply(k_means[,-1], 1, which.max) %>% c(k_range) %>% unique()
  
  dat <- dat[, k_order]
  colnames(dat) <- sprintf("K%02d", k_range)
  bind_cols(fam, dat)
}

read_log <- function(file) {
  lines <- read_lines(file)
  like <- str_split(lines[str_starts(lines, "Loglikelihood: ")], ":")[[1]][2]
  cver <- str_split(lines[str_starts(lines, "CV error")], ":")[[1]][2]
  tibble(
    likelihood = as.numeric(like),
    cv_error = as.numeric(cver),
  )
}

plot_distruct <- function(dat, group, ..., K = NULL) {
  if(!rlang::quo_is_null(enquo(K))) dat <- group_by(dat, {{K}})

  dat <- dat %>%
    arrange(..., {{group}}) %>%
    mutate(i = seq_len(n())) %>%
    ungroup() %>%
    pivot_longer(
      starts_with("K") & where(is.numeric),
      names_to = "cluster",
      values_to = "rate"
    ) %>%
    filter(!is.na(rate))
  
  dat_axis <- dat %>%
    group_by({{group}}, ...) %>%
    summarise(bar = max(i) + .5, pos = mean(i))
  
  ggplot(dat, aes(i, rate, fill = cluster)) +
    geom_col(width = 1) +
    geom_vline(
      aes(xintercept = bar), size = .1, color = "#454545",
      data = dat_axis
    ) +
    scale_x_continuous(
      "", expand = c(0, 0),
      breaks = pull(dat_axis, pos),
      labels = pull(dat_axis, {{group}})
    ) +
    scale_y_continuous(
      "", expand = c(0, 0),
      breaks = c(0, 1),
      labels = scales::percent
    ) +
    theme(
      axis.text.x = element_text(size = 6, angle = -90, hjust = 0, vjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(fill = NA),
      panel.spacing.x = unit(2, "pt"),
      panel.spacing.y = unit(2, "pt"),
      legend.position = "none"
    )
}