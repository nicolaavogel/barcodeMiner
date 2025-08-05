setwd("/Users/bfj994/Documents/barcodeMiner/")
# Install and load required package
library(VennDiagram)
library(grid)
library(rentrez)
library(dplyr)
library(ggtree)
library(ape)
library(taxize)
library(treeio)
library(tidyr)



# Load your accession lists
full_acc <- read.csv("full_accesion_list.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]
trnL_acc <- read.csv("trnL_accession_list.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]
barcode_acc <- read.csv("barcode_accession_list.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]

# Make the named list
venn_list <- list(
  Full = full_acc,
  trnL = trnL_acc,
  Barcode = barcode_acc
)

# Draw Venn Diagram
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("Full", "trnL", "Barcode"),
  filename = NULL,  # Use NULL to plot to R graphics device
  output = TRUE,
  fill = c("red", "green", "blue"),
  alpha = 0.5,
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.1
)

# Display it
grid.newpage()
grid.draw(venn.plot)

# Read names.dmp
names_df <- read.table("names.dmp", sep = "|", quote = "", stringsAsFactors = FALSE, fill = TRUE, header = FALSE)
names_df <- names_df[, c(1, 2, 4)]
colnames(names_df) <- c("tax_id", "name_txt", "name_class")
#names_df <- trimws(names_df)

# Keep only scientific names
names_df <- subset(names_df, name_class == "\tscientific name\t")

# Read nodes.dmp
nodes_df <- read.table("nodes.dmp", sep = "|", quote = "", stringsAsFactors = FALSE, fill = TRUE, header = FALSE)
nodes_df <- nodes_df[, c(1, 2, 3)]
colnames(nodes_df) <- c("tax_id", "parent_tax_id", "rank")
#nodes_df <- trimws(nodes_df)
acc_taxid <- read.table("refseq_plastid.genomic.acc2taxid.gz", header = TRUE, sep = "\t", quote = "")
# Filter only the accessions you care about

acc_taxid <- subset(acc_taxid, accession %in% full_acc)
get_lineage <- function(tax_id, nodes, names_df) {
  lineage <- list()
  while (TRUE) {
    node <- nodes[nodes$tax_id == tax_id, ]
    name <- names_df[names_df$tax_id == tax_id, "name_txt"]
    
    if (nrow(node) == 0 || length(name) == 0 || is.na(name[1])) break
    
    rank_clean <- trimws(node$rank)
    name_clean <- trimws(name[1])
    
    if (rank_clean != "") {
      lineage[[rank_clean]] <- name_clean
    }
    
    if (tax_id == node$parent_tax_id) break
    tax_id <- node$parent_tax_id
  }
  return(lineage)
}

accessions <- unique(acc_taxid$accession)
annotated <- data.frame()

for (i in seq_along(accessions)) {
  acc <- accessions[i]
  tax_id <- acc_taxid$taxid[acc_taxid$accession == acc]
  print(tax_id)
  lineage <- get_lineage(tax_id, nodes_df, names_df)
  print(lineage)
  annotated <- rbind(annotated, data.frame(
    accession = acc,
    species   = if (!is.null(lineage[["species"]])) lineage[["species"]] else NA,
    genus     = if (!is.null(lineage[["genus"]])) lineage[["genus"]] else NA,
    family    = if (!is.null(lineage[["family"]])) lineage[["family"]] else NA,
    order     = if (!is.null(lineage[["order"]])) lineage[["order"]] else NA,
    class     = if (!is.null(lineage[["class"]])) lineage[["class"]] else NA,
    phylum    = if (!is.null(lineage[["phylum"]])) lineage[["phylum"]] else NA,
    stringsAsFactors = FALSE
  ))
}

annotated$source <- NA

annotated$source[annotated$accession %in% full_acc] <- "full"
annotated$source[annotated$accession %in% trnL_acc & annotated$accession %in% full_acc] <- "full+trnL"
annotated$source[annotated$accession %in% barcode_acc & annotated$accession %in% full_acc] <- "full+barcode"
annotated$source[annotated$accession %in% full_acc &
                   annotated$accession %in% trnL_acc &
                   annotated$accession %in% barcode_acc] <- "all"

# Fill in remaining:
annotated$source[is.na(annotated$source)] <- "other"
annotated <- annotated %>%
  left_join(acc_taxid, by = "accession")
# Select taxids for which you want to build the tree (only those from your annotations)
taxid_set <- unique(annotated$taxid)

#write.csv(annotated, "new_full_annotated_save.tsv")
annotated <- read.csv("new_full_annotated_save.tsv", sep = ",", header =T)
# Create edge list tracing up to root
get_path_to_root <- function(taxid, nodes_df) {
  path <- c()
  while (!is.na(taxid) && taxid != 1) {
    parent <- nodes_df$parent_tax_id[nodes_df$tax_id == taxid]
    if (length(parent) == 0) break
    path <- rbind(path, data.frame(child = taxid, parent = parent))
    taxid <- parent
  }
  return(path)
}
############## ALL SPECIES #####################
all_edges <- do.call(rbind, lapply(taxid_set, get_path_to_root, nodes_df = nodes_df))
all_edges <- distinct(all_edges)

# Create tree from edge list
tree_tbl <- as_tibble(all_edges)
tree <- as.phylo(tree_tbl)

# Keep tree$tip.label as taxids (do NOT replace them)
tree$tip.label <- as.character(tree$tip.label)  # if not already

tip_data <- data.frame(taxid = as.integer(tree$tip.label)) %>%
  left_join(names_df[names_df$name_class == "\tscientific name\t", c("tax_id", "name_txt")],
            by = c("taxid" = "tax_id")) %>%
  left_join(annotated %>% select(taxid, source) %>% distinct(), by = "taxid")

# Plot with colored branches and tips
p <- ggtree(tree, aes(color = source)) %<+% tip_data +
  geom_tree(aes(color = source), size = 0.8) +        # color branches
  geom_tippoint(size = 2) +                           # color tips (inherits color)
  geom_tiplab(aes(label = name_txt), hjust = -0.1, size = 3) +
  theme_tree2() +
  labs(title = "Tree with Tips and Branches Colored by Source") +
  scale_color_brewer(palette = "Set1", na.value = "gray")

print(p)

############# EXCLUDING SPECIES THAT ARE PRESENT IN ALL DBs ####################

annotated_filtered <- annotated %>%
  filter(source != "all")
taxid_set1 <- unique(annotated_filtered$taxid)

# Get new edge list from filtered taxids
all_edges1 <- do.call(rbind, lapply(taxid_set1, get_path_to_root, nodes_df = nodes_df))
all_edges1 <- distinct(all_edges1)

# Create tree
tree_tbl1 <- as_tibble(all_edges1)
tree1 <- as.phylo(tree_tbl1)
tip_data1 <- data.frame(taxid = as.integer(tree1$tip.label)) %>%
  left_join(names_df[names_df$name_class == "\tscientific name\t", c("tax_id", "name_txt")],
            by = c("taxid" = "tax_id")) %>%
  left_join(annotated_filtered %>% select(taxid, source) %>% distinct(), by = "taxid")
ggtree(tree1, aes(color = source)) %<+% tip_data +
  geom_tree(aes(color = source), size = 0.8) +
  geom_tippoint(size = 2) +
  geom_tiplab(aes(label = name_txt), hjust = -0.1, size = 3) +
  theme_tree2() +
  labs(title = "Tree with Source-Colored Branches & Tips (Excluding 'all')") +
  scale_color_brewer(palette = "Set1", na.value = "gray")

###### full families only ######################################################
# Clean up
annotated2 <- annotated %>%
  mutate(source = trimws(source),
         species = trimws(species),
         family = trimws(family),
         genus = trimws(genus))

# Count unique sources per family
family_source_counts <- annotated2 %>%
  group_by(family) %>%
  summarise(
    n_sources = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ", ")
  )
# Count number of references per family
family_reference_counts <- annotated2 %>%
  group_by(family) %>%
  summarise(n_references = n(), .groups = "drop")


########################################################################
full_only_families <- family_source_counts %>%
  filter(n_sources == 1 & sources == "full") %>%
  pull(family)

trnL_only_families <- family_source_counts %>%
  filter(n_sources == 1 & sources == "full+trnL") %>%
  pull(family)

barcode_only_families <- family_source_counts %>%
  filter(n_sources == 1 & sources == "full+barcode") %>%
  pull(family)

# Combine all families into one vector with a source label
all_families_df <- tibble(
  family = c(full_only_families, trnL_only_families, barcode_only_families),
  source = c(
    rep("No trnL information or barcode region found", length(full_only_families)),
    rep("Not found as a barcode region, but trnL information present", length(trnL_only_families)),
    rep("Not found with trnL information, but found as a barcode region", length(barcode_only_families))
  )
)

# Get taxids for all families at once
family_taxids <- get_uid(all_families_df$family, ask = FALSE)

# Map taxid back to family for all
family_taxids_vec <- sapply(family_taxids, identity)
names(family_taxids_vec) <- names(family_taxids_vec)

taxid_to_family <- setNames(names(family_taxids_vec), as.character(family_taxids_vec))

# Get classification for all taxids
family_lineages <- classification(family_taxids_vec, db = "ncbi")

# Add family and taxid info to each lineage df
family_lineages_named <- lapply(names(family_lineages), function(tid) {
  df <- family_lineages[[tid]]
  df$family <- taxid_to_family[tid]
  df$taxid <- tid
  df
})

lineage_filtered <- bind_rows(family_lineages_named) %>%
  filter(rank %in% c("phylum", "family")) %>%
  select(taxid, rank, name, family) %>%
  rename(family_name = family) %>%  # rename to avoid conflict
  distinct() %>%
  pivot_wider(names_from = rank, values_from = name) %>%
  relocate(taxid, phylum, family, family_name)

# Join with all_families_df to get the source column for each family
lineage_filtered <- lineage_filtered %>%
  left_join(all_families_df, by = "family")
lineage_filtered <- lineage_filtered %>%
  left_join(family_reference_counts, by = "family")

ggplot(lineage_filtered, aes(x = phylum, y = family, color = source, size = n_references)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(2, 10)) +  # Adjust dot size range as needed
  theme_bw(base_size = 15) +
  labs(
    color = "Database",
    size = "Number of References",
    x = "Phylum",
    y = "Family"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Family_overview_in_DBs.png", height = 23, width = 17)

##############full genera only####################################

genus_source_counts <- annotated2 %>%
  group_by(genus) %>%
  summarise(
    n_sources = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ", ")
  )

full_only_genera <- genus_source_counts %>%
  filter(n_sources == 1 & sources == "full")

print(full_only_genera)

####################################################################

full_lin <- bind_rows(family_lineages_named) %>%
  filter(rank %in% c("phylum", "order", "class", "family")) %>%
  select(taxid, rank, name, family) %>%
  rename(family_name = family) %>%  # rename to avoid conflict
  distinct() %>%
  pivot_wider(names_from = rank, values_from = name) %>%
  relocate(taxid, phylum, family, family_name)

# Join with all_families_df to get the source column for each family
full_lin <- full_lin %>%
  left_join(all_families_df, by = "family")
full_lin <- full_lin %>%
  left_join(family_reference_counts, by = "family")

full_lin_filtered <- full_lin |>
  filter(phylum == "Streptophyta")

ggplot(full_lin_filtered, aes(x = class, y = family, color = source, size = n_references)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(2, 10)) +  # Adjust dot size range as needed
  theme_bw(base_size = 15) +
  labs(
    color = "Database",
    size = "Number of References",
    x = "Class",
    y = "Family"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Family_overview_streptophyta_class.png", height = 23, width = 17)


################################################################################
#### new visuals 
# First gather the data by taxonomic rank and source
df_long <- annotated |>
  pivot_longer(cols = c(family, order, class, phylum), names_to = "rank", values_to = "rank_name")

# Count combinations
rank_source_summary <- df_long |>
  group_by(rank, rank_name, source) |>
  summarise(n = n(), .groups = "drop")

# Plot: For each rank, how often each source appears
ggplot(rank_source_summary, aes(x = reorder(rank_name, -n), y = n, fill = source)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_log10() +
  facet_wrap(~rank, scales = "free_x", ncol = 1) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Taxon Name", y = "Count", fill = "Source Type",
       title = "Distribution of Source Types Across Taxonomic Ranks")
ggsave("DB_distribution.png", height = 10, width = 20)

####
# Load required libraries
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# Your data: assuming it's already in a dataframe called `annotated`
# Step 1: Convert taxonomic levels to numeric features for clustering
taxonomy_levels <- c("phylum", "class", "order", "family", "genus", "species")

# Create a numeric matrix based on taxonomy (each column = one-hot encoding of taxonomic label)
taxonomy_matrix <- annotated %>%
  select(all_of(taxonomy_levels)) %>%
  mutate(across(everything(), as.factor)) %>%
  model.matrix(~ . - 1, data = .)  # one-hot encoding

# Step 2: Generate row annotation for source
# Assign each source a color
source_colors <- RColorBrewer::brewer.pal(length(unique(annotated$source)), "Set3")
names(source_colors) <- unique(annotated$source)

annotation_row <- data.frame(source = annotated$source)
rownames(annotation_row) <- annotated$accession
ann_colors <- list(source = source_colors)

# Step 3: Create heatmap
pheatmap(
  mat = taxonomy_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 10,
  main = "Phylogenetic Clustering by Taxonomy\nColor-coded by Source"
)
