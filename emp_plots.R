setwd("/Users/bfj994/Documents/barcodeMiner/emp_results/")
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
library(purrr)

dir_path <- "/Users/bfj994/Documents/barcodeMiner/emp_results/"
lca_files <- list.files(path = dir_path, pattern = "\\.lca\\.gz$", full.names = TRUE)
# Read and combine all files with an extra column indicating the base filename
emp_lca_data <- rbindlist(lapply(lca_files, function(file) {
  df <- fread(file)  # use read.table(file, header = TRUE) if not using data.table
  file_label <- sub("\\.lca\\.gz$", "", basename(file))
  df$sample <- file_label
  return(df)
}))
emp_lca_data <- as.data.table(emp_lca_data)
emp_lca_data[, c("taxid", "name", "rank") := tstrsplit(lca, ':"|":"', perl=TRUE)]

emp_lca <- emp_lca_data |>
  mutate(
    sample_clean = str_remove(sample, "\\.lca\\.gz$"),
    parts = str_split(sample_clean, "_"),
    sample = map_chr(parts, 1),
    db = map_chr(parts, 4),
    rank_name = paste0(rank, ": ", name)
  ) |>
  select(-parts, -sample_clean)

emp_lca |>
  filter(db == "barcode") |>
  filter(len > 19)|>
  group_by(sample, rank_name) |>
  summarise(n = n(), .groups = "drop") |>
  filter(!is.na(rank_name)) |>
  ggplot(aes(x = reorder(rank_name, -n), y = n)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~ sample, scales = "free_y") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "LCA Assignment (Rank: Name)",
    y = "Number of Reads",
    title = "LCA Assignments for 'barcode' DB"
  )
ggsave("emp_data_barcode.png")

emp_lca |>
  filter(db == "trnl") |>
  filter(len > 19)|>
  group_by(sample, rank_name) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n >= 5) |>
  filter(rank_name != "species\": Andrographis paniculata") |>
  ggplot(aes(x = reorder(rank_name, -n), y = n)) +
  geom_col(fill = "darkgreen") +
  facet_wrap(~ sample, scales = "free_y") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "LCA Assignment (Rank: Name)",
    y = "Number of Reads",
    title = "LCA Assignments for 'trnL' DB"
  )
ggsave("emp_data_trnL.png")
