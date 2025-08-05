setwd("/Users/bfj994/Documents/barcodeMiner/emp_results/")
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
library(purrr)
library(tidyr)


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
  # Filter out repeated header rows
  filter(
    queryid != "queryid",
    lca != "lca",
    taxa_path != "taxa_path"
  ) |>
  # Clean and parse the sample column
  mutate(
    sample_clean = str_remove(sample, "\\.lca\\.gz$"),
    parts = str_split(sample_clean, "_"),
    sample = map_chr(parts, 1),         # e.g., "CGG3016617" or "Lib"
    db = map_chr(parts, ~ tail(., 1)),  # take last part, e.g., "barcode"
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
  facet_wrap(~ sample, scales = "free_y", ncol = 4) +
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
ggsave("emp_data_barcode.png", width = 20, height = 20)

emp_lca |>
  filter(db == "trnl") |>
  filter(len > 19)|>
  group_by(sample, rank_name) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n >= 5) |>
  filter(rank_name != "species\": Andrographis paniculata") |>
  ggplot(aes(x = reorder(rank_name, -n), y = n)) +
  geom_col(fill = "darkgreen") +
  facet_wrap(~ sample, scales = "free_y", ncol = 4) +
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
ggsave("emp_data_trnL.png",width = 20, height = 20)

emp_lca |>
  filter(db == "full") |>
  filter(len > 29)|>
  group_by(sample, rank_name) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n >= 500) |>
  filter(rank_name != "species\": Andrographis paniculata") |>
  ggplot(aes(x = reorder(rank_name, -n), y = n)) +
  geom_col(fill = "darkred") +
  facet_wrap(~ sample, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "LCA Assignment (Rank: Name)",
    y = "Number of Reads",
    title = "LCA Assignments for 'full' DB"
  )
ggsave("emp_data_full.png", width = 20, height = 20)


################################################################################
###### Difference 

library(tidyr)
library(dplyr)
library(ggplot2)

# Step 1: Combine counts as before
summarise_db <- function(df, db_name, min_len, min_n = 1) {
  df |>
    filter(db == db_name, len > min_len) |>
    group_by(sample, rank_name) |>
    summarise(n = n(), .groups = "drop") |>
    filter(n >= min_n) |>
    mutate(db = db_name)
}

barcode_counts <- summarise_db(emp_lca, "barcode", 19)
trnl_counts    <- summarise_db(emp_lca, "trnl", 19, min_n = 5)
full_counts    <- summarise_db(emp_lca, "full", 29, min_n = 1000)

combined_counts <- bind_rows(barcode_counts, trnl_counts, full_counts)

# Step 2: Pivot to wide format per sample
taxa_wide <- combined_counts |>
  pivot_wider(
    names_from = db, values_from = n, values_fill = 0
  )
# Add column counting number of DBs with count > 0
taxa_with_2plus <- taxa_wide |>
  mutate(
    db_nonzero = (barcode > 0) + (trnl > 0) + (full > 0)
  ) |>
  filter(db_nonzero >= 2) |>
  select(-db_nonzero)

# Convert to long format for plotting
plot_data_2plus <- taxa_with_2plus |>
  pivot_longer(cols = c(barcode, trnl, full),
               names_to = "db", values_to = "count") |>
  filter(count > 0)

# Plot
ggplot(plot_data_2plus, aes(x = reorder(rank_name, -count), y = count, fill = db)) +
  geom_col(position = "dodge") +
  facet_wrap(~ sample, scales = "free_y", ncol = 4) +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Taxon (Rank: Name)",
    y = "Read Count (log10 scale)"
  )
ggsave("taxa_shared_log10.png", width = 21, height = 18)

taxa_only_full <- taxa_wide |>
  filter(barcode == 0, trnl == 0, full > 0)

plot_data_only_full <- taxa_only_full |>
  select(sample, rank_name, full) |>
  rename(count = full) |>
  mutate(db = "full")

ggplot(plot_data_only_full, aes(x = reorder(rank_name, -count), y = count)) +
  geom_col(fill = "darkred") +
  facet_wrap(~ sample, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Taxon (Rank: Name)",
    y = "Read Count",
    title = "Taxa Only Found in 'full' DB per Sample"
  )
ggsave("taxa_only_full.png", width = 21, height = 20)

################################################################################
emp_lca <- emp_lca_data |>
  # Filter out repeated header rows
  filter(
    queryid != "queryid",
    lca != "lca",
    taxa_path != "taxa_path"
  ) |>
  # Clean and parse the sample column
  mutate(
    sample_clean = str_remove(sample, "\\.lca\\.gz$"),
    parts = str_split(sample_clean, "_"),
    sample = map_chr(parts, 1),         # e.g., "CGG3016617" or "Lib"
    db = map_chr(parts, ~ tail(., 1)),  # take last part, e.g., "barcode"
    rank_name = paste0(rank, ": ", name)
  ) |>
  select(-parts, -sample_clean)

###############################################################################

extract_taxonomy <- function(lineage_string) {
  # Split the string and extract names and ranks using regex
  entries <- strsplit(lineage_string, ";")[[1]]
  parsed <- stringr::str_match(entries, '^[^:]+:"([^"]+)":"([^"]+)"$')
  
  # Clean up parsing result
  tax_names <- parsed[, 2]
  tax_ranks <- tolower(parsed[, 3])  # normalize ranks to lowercase
  
  # Filter out NAs from parsing errors
  valid <- !is.na(tax_names) & !is.na(tax_ranks)
  tax_names <- tax_names[valid]
  tax_ranks <- tax_ranks[valid]
  
  # Make a named vector
  taxonomy <- setNames(tax_names, tax_ranks)
  
  # Use safe access (check if the rank exists first)
  get_rank <- function(rank) {
    if (rank %in% names(taxonomy)) taxonomy[[rank]] else NA_character_
  }
  
  tibble(
    genus  = get_rank("genus"),
    family = get_rank("family"),
    order  = get_rank("order"),
    class  = get_rank("class"),
    phylum = get_rank("phylum")
  )
}




############################################################################

emp_lca2 <- emp_lca |>
  mutate(
    genus  = if_else(rank == 'genus"', name, NA_character_),
    family = if_else(rank == 'family"', name, NA_character_),
    class  = if_else(rank == 'class"', name, NA_character_)
  )
summarise_taxa_rank <- function(df, tax_rank, min_len, min_n = 1) {
  df |>
    filter(len > min_len) |>
    filter(!is.na(.data[[tax_rank]])) |>
    group_by(sample, db, !!sym(tax_rank)) |>
    summarise(n = n(), .groups = "drop") |>
    filter(n >= min_n) |>
    rename(taxon = !!sym(tax_rank))
}
genus_counts  <- summarise_taxa_rank(emp_lca2, "genus", 19)
family_counts <- summarise_taxa_rank(emp_lca2, "family", 19)
class_counts  <- summarise_taxa_rank(emp_lca2, "class", 19)

compare_taxa_across_dbs <- function(taxa_df, tax_rank) {
  taxa_df_wide <- taxa_df |>
    pivot_wider(
      names_from = db,
      values_from = n,
      values_fill = 0
    )
  
  # Count how many DBs found this taxon
  taxa_df_wide <- taxa_df_wide |>
    mutate(
      db_nonzero = rowSums(across(c(barcode, trnl, full), ~ . > 0))
    )
  
  # Plot
  taxa_df_long <- taxa_df_wide |>
    pivot_longer(cols = c(barcode, trnl, full),
                 names_to = "db", values_to = "count") |>
    filter(count > 0)
  
  p <- ggplot(taxa_df_long, aes(x = reorder(taxon, -count), y = count, fill = db)) +
    geom_col(position = "dodge") +
    facet_wrap(~ sample, scales = "free_y", ncol = 4) +
    scale_y_log10() +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.text = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = paste("Taxon (", tax_rank, ")", sep = ""),
      y = "Read Count (log10 scale)",
      title = paste("Taxa (", tax_rank, ") Found Across DBs", sep = "")
    )
  
  ggsave(paste0("taxa_", tax_rank, "_across_dbs.png"), p, width = 20, height = 18)
}

compare_taxa_across_dbs(genus_counts, "genus")
ggsave("genus_comp_all_DBs.png", height = 20, width = 20)
compare_taxa_across_dbs(family_counts, "family")
ggsave("family_comp_all_DBs.png", height = 20, width = 20)






####### per sample plots GENUS
genus_counts <- emp_lca |>
  filter(db != "ex") |>
  filter(rank == 'genus"', !is.na(name)) |>
  count(sample, db, name)

# Pivot to wide format
genus_wide <- genus_counts |>
  pivot_wider(
    names_from = db,
    values_from = n,
    values_fill = 0
  )

# Function to plot per sample
plot_sample_taxa <- function(df, sample_id) {
  df_sample <- df |> filter(sample == sample_id)
  
  df_long <- df_sample |>
    pivot_longer(
      cols = -c(sample, name),
      names_to = "db",
      values_to = "count"
    )
  df_long <- df_long |> filter(count > 0)
  
  ggplot(df_long, aes(x = name, y = count, fill = db)) +
    geom_col(position = "dodge") +
    scale_y_log10()+
    labs(
      title = paste("Genus composition for", sample_id),
      x = "Genus",
      y = "Read count (log10)",
      fill = "Database"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1))
}

# Create plot directory if it doesn't exist
if (!dir.exists("plots")) dir.create("plots")

# Generate plots
sample_ids <- unique(genus_wide$sample)
sample_plots <- map(sample_ids, ~ plot_sample_taxa(genus_wide, .x))
names(sample_plots) <- sample_ids

# Save plots to disk
walk2(sample_ids, sample_plots, ~ ggsave(
  filename = paste0("plots/", .x, "_genus_plot.png"),
  plot = .y,
  width = 20, height = 15
))

# Optionally, preview the first plot
sample_plots[[7]]

########### Per sample plots FAMILY
####### per sample plots 
genus_counts <- emp_lca |>
  filter(db != "ex") |>
  filter(rank == 'family"', !is.na(name)) |>
  count(sample, db, name)

# Pivot to wide format
genus_wide <- genus_counts |>
  pivot_wider(
    names_from = db,
    values_from = n,
    values_fill = 0
  )

# Function to plot per sample
plot_sample_taxa <- function(df, sample_id) {
  df_sample <- df |> filter(sample == sample_id)
  
  df_long <- df_sample |>
    pivot_longer(
      cols = -c(sample, name),
      names_to = "db",
      values_to = "count"
    )
  df_long <- df_long |> filter(count > 0)
  
  ggplot(df_long, aes(x = name, y = count, fill = db)) +
    geom_col(position = "dodge") +
    scale_y_log10()+
    labs(
      title = paste("Family composition for", sample_id),
      x = "Family",
      y = "Read count (log10)",
      fill = "Database"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1))
}

# Create plot directory if it doesn't exist
if (!dir.exists("plots")) dir.create("plots")

# Generate plots
sample_ids <- unique(genus_wide$sample)
sample_plots <- map(sample_ids, ~ plot_sample_taxa(genus_wide, .x))
names(sample_plots) <- sample_ids

# Save plots to disk
walk2(sample_ids, sample_plots, ~ ggsave(
  filename = paste0("plots/", .x, "_family_plot.png"),
  plot = .y,
  width = 20, height = 15
))

# Optionally, preview the first plot
sample_plots[[4]]

###############################################################################
#full data 
library(data.table)

# Convert to data.table for speed
setDT(emp_lca)

# Use vectorized string splitting to extract parts
# This will be WAY faster if you only care about known ranks

extract_ranks_fast <- function(x, ranks = c("genus", "family", "order", "class", "phylum")) {
  out <- vector("list", length(ranks))
  names(out) <- ranks
  
  for (rank in ranks) {
    pattern <- paste0(':"([^"]+)":"', rank, '"')
    match <- stringr::str_match(x, pattern)[, 2]
    out[[rank]] <- match
  }
  
  as.data.table(out)
}

# Apply to full column (this is vectorized!)
tax_split <- extract_ranks_fast(emp_lca$taxa_path)

# Bind the result back to the original table
emp_lca_with_tax <- cbind(emp_lca, tax_split)

######## SINGLE SAMPLE 

test <- emp_lca |>
  filter(sample == "CGG3017472") |>
  mutate(tax_split = map(taxa_path, extract_taxonomy)) |> 
  unnest(tax_split)

library(dplyr)
library(ggplot2)
library(forcats)


genus_counts <- test %>%
  filter(db != "ex") %>%
  filter(!is.na(genus)) %>%
  count(genus, db, sort = TRUE) %>%
  filter(n > 0) %>%
  group_by(genus) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(
    log_n = log10(n),
    genus = fct_reorder(genus, total)
  )

# Family counts
family_counts <- test %>%
  filter(db != "ex") %>%
  filter(!is.na(family)) %>%
  count(family, db, sort = TRUE) %>%
  filter(n > 0) %>%
  group_by(family) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(
    log_n = log10(n),
    family = fct_reorder(family, total)
  )
ggplot(genus_counts, aes(x = genus, y = log_n, fill = db)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Log10 Counts of Genera by Database for sample CGG3017472",
    x = "Genus",
    y = "log10(Count)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


ggplot(family_counts, aes(x = family, y = log_n, fill = db)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Log10 Counts of Families by Database",
    x = "Family",
    y = "log10(Count)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Count class-level occurrences
class_counts <- test %>%
  filter(db != "ex") %>%
  filter(!is.na(class)) %>%
  count(class, db, sort = TRUE) %>%
  filter(n > 0) %>%
  group_by(class) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(
    log_n = log10(n),
    class = fct_reorder(class, total)
  )

# Plot
ggplot(class_counts, aes(x = class, y = log_n, fill = db)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Log10 Counts of Classes by Database for sample CGG3017472",
    x = "Class",
    y = "log10(Count)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")


# Find top 30 families
top_families <- test %>%
  filter(db != "ex") %>%
  filter(!is.na(family)) %>%
  count(family, sort = TRUE) %>%
  slice_max(n, n = 30) %>%
  pull(family)

# Count by family and db, only for top families
family_top_counts <- test %>%
  filter(db != "ex") %>%
  filter(family %in% top_families, !is.na(db)) %>%
  count(family, db) %>%
  filter(n > 0) %>%
  mutate(
    family = fct_reorder(family, n)
  )
ggplot(family_top_counts, aes(x = family, y = n, fill = db)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~db, scales = "free") +
  labs(
    title = "Read Counts of Top 30 Families by Database for sample CGG3017472",
    x = "Family",
    y = "Read Count"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )
ggsave("CGG3017472_top_30_fam.png")

# Find top 30 families
top_genera <- test %>%
  filter(db != "ex") %>%
  filter(!is.na(genus)) %>%
  count(genus, sort = TRUE) %>%
  slice_max(n, n = 30) %>%
  pull(genus)

# Count by family and db, only for top families
genera_top_counts <- test %>%
  filter(db != "ex") %>%
  filter(genus %in% top_genera, !is.na(db)) %>%
  count(genus, db) %>%
  filter(n > 0) %>%
  mutate(
    genus = fct_reorder(genus, n)
  )
ggplot(genera_top_counts, aes(x = genus, y = n, fill = db)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~db, scales = "free") +
  labs(
    title = "Read Counts of Top 30 Genera by Database for sample CGG3017472",
    x = "Genus",
    y = "Read Count"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )
ggsave("CGG3017472_top_30_gen.png")
