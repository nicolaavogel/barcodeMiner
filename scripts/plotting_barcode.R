setwd("/Users/bfj994/Documents/barcodeMiner/results/")

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
library(purrr)

dir_path <- "/Users/bfj994/Documents/barcodeMiner/results/"
lca_files <- list.files(path = dir_path, pattern = "\\.lca\\.gz$", full.names = TRUE)
# Read and combine all files with an extra column indicating the base filename
lca_data <- rbindlist(lapply(lca_files, function(file) {
  df <- fread(file)  # use read.table(file, header = TRUE) if not using data.table
  file_label <- sub("\\.lca\\.gz$", "", basename(file))
  df$sample <- file_label
  return(df)
}))
lca_data <- as.data.table(lca_data)

# Split LCA into taxid, name, and rank
lca_data[, c("taxid", "name", "rank") := tstrsplit(lca, ':"|":"', perl=TRUE)]


lca <- lca_data |>
  mutate(
    # Remove file extension if present
    sample_clean = str_remove(sample, "\\.lca\\.gz$"),
    
    # Split the sample into components
    parts = str_split(sample_clean, "_"),
    
    # Extract amount (2nd element)
    amount = case_when(
      str_detect(map_chr(parts, 2), "k$") ~ as.numeric(str_remove(map_chr(parts, 2), "k")) * 1e3,
      str_detect(map_chr(parts, 2), "m$") ~ as.numeric(str_remove(map_chr(parts, 2), "m")) * 1e6,
      TRUE ~ NA_real_
    ),
    
    # Compute coverage
    coverage = (amount * 45) / 160000,
    
    # Extract bootstrap (last element)
    bootstrap = as.integer(map_chr(parts, ~ .x[length(.x)])),
    
    # Extract db (everything between amount and bootstrap)
    db = map_chr(parts, ~ if (length(.x) > 3) paste(.x[3:(length(.x) - 1)], collapse = "_") else .x[3]),
    
    # Create rank label
    rank_name = paste0(rank, ": ", name)
  ) |>
  select(-parts, -sample_clean)


lca <- lca |>
  filter(!is.na(name) & !is.na(rank)) |>
  mutate(
    rank_name = paste0(rank, ": ", name)
  )

# Step 2: Simulated total reads per coverage
coverage_reads_lookup <- lca |>
  distinct(coverage, amount) |>
  deframe()  # turns into named vector: coverage -> amount

# Step 3: Consistent color mapping for DBs
db_levels <- unique(lca$db)
db_colors <- setNames(RColorBrewer::brewer.pal(n = length(db_levels), name = "Set2"), db_levels)

rank_summary <- lca |>
  filter(db != "trnl_ex") |>
  filter(len > 19) |>
  group_by(coverage, rank_name, db, bootstrap) |>
  summarise(n = n(), .groups = "drop") |>  # this gives "n per simulation"
  group_by(coverage, rank_name, db) |>
  summarise(
    mean_n = mean(n),
    se_n = sd(n) / sqrt(n()),
    ci75 = 0.674 * se_n,
    ymin = pmax(mean_n - ci75, 1e-3),
    ymax = mean_n + ci75,
    .groups = "drop"
  ) |>
  filter(mean_n > 10)

# Plotting with log10 scale, no zero points included
plots_by_coverage <- rank_summary |>
  group_split(coverage) |>
  map(function(df) {
    coverage_val <- unique(df$coverage)
    sim_reads <- coverage_reads_lookup[as.character(coverage_val)]
    
    ggplot(df, aes(x = reorder(rank_name, -mean_n), y = mean_n, color = db)) +
      geom_point(position = position_dodge(width = 0.7), size = 3, alpha = 0.9) +
      scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        expand = expansion(mult = c(0.01, 0.1))
      ) +
      scale_color_manual(values = db_colors) +
      theme_bw(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.spacing = unit(1, "lines")
      ) +
      labs(
        title = paste0(round(coverage_val), "x coverage (Simulated reads: ", format(sim_reads, big.mark = ","), ")"),
        x = "Taxonomic Rank: Name",
        y = "Mean Number of Assignments (log10 scale)",
        color = "Method (DB)"
      )
    
  })

# Save plots
walk2(plots_by_coverage, seq_along(plots_by_coverage), function(plot, i) {
  ggsave(
    filename = paste0("BN_rank_assignments_", i, ".png"),
    plot = plot,
    width = 12, height = 8, dpi = 300
  )
})

