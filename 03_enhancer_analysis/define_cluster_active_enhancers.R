library(dplyr)
library(readr)

count_dir <- "results/cluster_enhancer_counts"
output_dir <- "results/cluster_active_enhancers"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  count_dir,
  pattern = "_enhancer_counts.txt$",
  full.names = TRUE
)

for (file in files) {

  cluster <- sub(
    "_enhancer_counts.txt$",
    "",
    basename(file)
  )

  df <- read_tsv(
    file,
    col_names = FALSE,
    show_col_types = FALSE
  )

  colnames(df)[1:3] <- c("chr", "start", "end")
  colnames(df)[ncol(df)] <- "count"

  df <- df %>%
    mutate(
      CPM = count / sum(count, na.rm = TRUE) * 1e6
    )

  active <- df %>%
    filter(CPM > 1) %>%
    select(chr, start, end)

  write_tsv(
    active,
    file.path(
      output_dir,
      paste0(cluster, "_active_enhancers_CPM1.bed")
    ),
    col_names = FALSE
  )
}
