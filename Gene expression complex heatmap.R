# Libraries
library(readxl)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Read the data
data <- read_excel("/Users/macbook/Desktop/New GW/Transcriptome Data/expression analysis.xlsx")

# Check for required columns
required_columns <- c("gene_id", "SMPCT-1", "SMPCT-2", "SMPCT-3", "GC", "cptmt_chr", "Type")
missing_columns <- setdiff(required_columns, colnames(data))
if (length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
}

# Remove any rows with NA values
data <- na.omit(data)

# Handle duplicate gene_ids
data <- data %>%
  group_by(gene_id) %>%
  slice(1) %>%
  ungroup()

# Prepare the main expression matrix
mat <- data %>%
  select(gene_id, `SMPCT-1`, `SMPCT-2`, `SMPCT-3`) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# Prepare GC content data
gc_data <- data %>%
  select(gene_id, GC) %>%
  column_to_rownames("gene_id")

# Prepare chromosome data
chr_data <- data %>%
  select(gene_id, cptmt_chr) %>%
  column_to_rownames("gene_id")

# Prepare Type data
type_data <- data %>%
  select(gene_id, Type) %>%
  column_to_rownames("gene_id")

# Ensure all datasets have the same rows in the same order
common_genes <- Reduce(intersect, list(rownames(mat), rownames(gc_data), rownames(chr_data), rownames(type_data)))
mat <- mat[common_genes, ]
gc_data <- gc_data[common_genes, , drop = FALSE]
chr_data <- chr_data[common_genes, , drop = FALSE]
type_data <- type_data[common_genes, , drop = FALSE]

# Define new color schemes
expr_colors <- colorRamp2(c(min(mat), mean(mat), max(mat)), c("#3B0F6F", "#56C667", "#FDE725"))
gc_colors <- colorRamp2(c(min(gc_data$GC), max(gc_data$GC)), c("#3B0F6F", "#FDE725"))
# Define colors for "Type"
type_colors <- structure(c("#32127A", "gold"), names = unique(type_data$Type))
# Define colors for chromosomes
chr_levels <- unique(chr_data$cptmt_chr)
chr_colors <- setNames(rainbow(length(chr_levels)), chr_levels)

# Convert chromosome data to a vector
chr_vector <- as.vector(chr_data$cptmt_chr)

# Create the heatmap with gene IDs and Type variable displayed
ht_list <- Heatmap(mat, name = "Expression",
                   col = expr_colors,
                   show_row_names = TRUE,
                   row_names_gp = gpar(fontsize = 8),
                   row_names_side = "left",
                   width = unit(8, "cm"),
                   column_title = "Genes Expression",
                   cluster_rows = TRUE,
                   cluster_columns = FALSE) +
  
  Heatmap(gc_data$GC, name = "GC content",
          col = gc_colors,
          width = unit(5, "mm")) +
  
  Heatmap(chr_vector, name = "Chromosome",
          col = chr_colors,
          width = unit(5, "mm")) +
  
  Heatmap(type_data$Type, name = "Type",
          col = type_colors,
          width = unit(5, "mm"))

# Draw the heatmap
draw(ht_list)

# Print summary information
cat("Number of genes:", nrow(mat), "\n")
cat("Number of unique chromosomes:", length(chr_levels), "\n")
cat("Chromosomes:", paste(chr_levels, collapse = ", "), "\n")
cat("GC content range:", min(gc_data$GC), "to", max(gc_data$GC), "\n")
cat("Unique Types:", paste(unique(type_data$Type), collapse = ", "), "\n")
# Save as PDF
pdf("heatmap_output.pdf", width = 10, height = 8)  # Adjust size if needed

# Draw the heatmap
draw(ht_list)

# Close the PDF device
dev.off()

