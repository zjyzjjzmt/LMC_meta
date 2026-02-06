#!/usr/bin/env Rscript

# ==============================================================================
# Script Name: 06_visualization.R
# Description: Generates figures for the manuscript including Maps, Upset plots,
#              PCoA analyses, and Heatmaps.
# Author: [Your Name]
# Date: [Current Date]
# Dependencies: ggplot2, ggmap, sp, maps, UpSetR, tidyverse, vegan, plyr, 
#               pairwiseAdonis, pheatmap
# ==============================================================================

# --- Setup Environment ---
# Function to load packages or install if missing
load_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# Load required packages
packages <- c("ggplot2", "ggmap", "sp", "maps", "UpSetR", "tidyverse", 
              "grid", "vegan", "plyr", "pheatmap")
invisible(lapply(packages, load_package))

# Install pairwiseAdonis from GitHub if not present
if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
  library(pairwiseAdonis)
}

# Define Input/Output Directories (Relative paths for portability)
INPUT_DIR <- "./07_Visualization_Data"
OUTPUT_DIR <- "./08_Figures"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Global Map Visualization
# ==============================================================================
message("Generating Figure 1: Global Map...")

tryCatch({
  # Load Data
  map_data_path <- file.path(INPUT_DIR, "data_collection_final.csv")
  if (!file.exists(map_data_path)) stop("Map data file not found!")
  
  mydata <- read.csv(map_data_path, header = TRUE)
  
  # Prepare Coordinates
  visit.x <- mydata$longtitude
  visit.y <- mydata$latitude
  
  # Base Map
  mapworld <- borders("world", colour = "gray50", fill = "white")
  
  # Plot
  mp <- ggplot() + 
    mapworld + 
    ylim(-60, 90) +
    geom_point(aes(x = visit.x, y = visit.y, size = mydata$number), 
               color = "darkorange", alpha = 0.8) +
    scale_size(range = c(2.5, 6)) + # Adjusted range for better visibility
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  ggsave(file.path(OUTPUT_DIR, "Figure1_Map.pdf"), plot = mp, width = 12, height = 6)
  message("  -> Saved Figure1_Map.pdf")
}, error = function(e) { message("  Error in Map plotting: ", e$message) })

# ==============================================================================
# 2. Upset Plot (ARG Intersections)
# ==============================================================================
message("Generating Figure 2: Upset Plot...")

tryCatch({
  upset_data_path <- file.path(INPUT_DIR, "all-upset-4.csv")
  if (!file.exists(upset_data_path)) stop("Upset data file not found!")
  
  data_table <- read.csv(upset_data_path, row.names = 1)
  
  # Ensure binary format (0/1)
  # data_table[] <- lapply(data_table, function(x) ifelse(x != 0, 1, 0))
  
  pdf(file.path(OUTPUT_DIR, "Figure2_Upset.pdf"), width = 10, height = 6, onefile = FALSE)
  upset(data_table, 
        sets = c("Other", "CAM", "CHM", "PM"), 
        keep.order = TRUE, 
        mainbar.y.label = "Intersection Size", 
        sets.x.label = "Set Size", 
        main.bar.color = "black", 
        matrix.color = "black", 
        sets.bar.color = "#E41A1C", # Red
        order.by = "freq", 
        shade.color = "gray88", 
        scale.intersections = "identity",
        text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.3)
  )
  dev.off()
  message("  -> Saved Figure2_Upset.pdf")
}, error = function(e) { message("  Error in Upset plotting: ", e$message) })

# ==============================================================================
# 3. PCoA Analysis and Plotting
# ==============================================================================
message("Generating Figure 3: PCoA Analysis...")

tryCatch({
  otu_path <- file.path(INPUT_DIR, "ARG-phase.txt")
  group_path <- file.path(INPUT_DIR, "group-phase.txt")
  
  if (!file.exists(otu_path) || !file.exists(group_path)) stop("PCoA input files not found!")
  
  # Load OTU table
  otu <- read.delim(otu_path, row.names = 1, check.names = FALSE)
  otu <- data.frame(t(otu))
  otu <- otu[rowSums(otu) > 0, ] # Remove empty rows
  
  # Load Group info
  group <- read.delim(group_path, stringsAsFactors = FALSE)
  
  # Calculate Distance (Bray-Curtis)
  distance <- vegdist(otu, method = 'bray')
  
  # PCoA
  pcoa <- cmdscale(distance, k = (nrow(otu) - 1), eig = TRUE)
  
  # Explanation percentages
  pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
  
  # Prepare Plot Data
  sample_site <- data.frame({pcoa$point})[1:2]
  names(sample_site) <- c('PCoA1', 'PCoA2')
  sample_site$names <- rownames(sample_site)
  sample_site <- merge(sample_site, group, by = 'names', all.x = TRUE)
  
  # Order factors
  sample_site$type <- factor(sample_site$type, 
                             levels = c('Feedstock', 'Heating', 'Thermophilic', 'Cooling', 'Maturation'))
  
  # Calculate Hull for polygons
  group_border <- ddply(sample_site, 'type', function(df) df[chull(df[[2]], df[[3]]), ])
  
  # Plot
  pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2, group = type)) +
    geom_vline(xintercept = 0, color = 'gray', linewidth = 0.4) + 
    geom_hline(yintercept = 0, color = 'gray', linewidth = 0.4) +
    geom_polygon(data = group_border, aes(fill = type), alpha = 0.3) +
    geom_point(aes(color = type), size = 3, alpha = 0.8) +
    scale_color_manual(values = c('#4E7FAE', '#7B9A53', '#A97B42', '#9E5162', '#9A7399')) +
    scale_fill_manual(values = c('#A3C1DD', '#C1D4A0', '#DDBF9C', '#D9A3AD', '#D0BDD0')) +
    labs(x = paste0('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), 
         y = paste0('PCoA2: ', round(100 * pcoa_eig[2], 2), '%')) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "right")
  
  ggsave(file.path(OUTPUT_DIR, "Figure3_PCoA.pdf"), plot = pcoa_plot, width = 8, height = 6)
  message("  -> Saved Figure3_PCoA.pdf")
  
  # --- Statistics (PERMANOVA & ANOVA) ---
  message("  Running Statistical Tests...")
  stats_file <- file.path(OUTPUT_DIR, "PCoA_Statistics.txt")
  sink(stats_file)
  
  cat("=== PERMANOVA (adonis2, Bray-Curtis) ===\n")
  dist_df <- sample_site
  rownames(dist_df) <- dist_df$names
  # Match rows between distance matrix and metadata
  valid_samples <- intersect(rownames(as.matrix(distance)), rownames(dist_df))
  dist_subset <- as.matrix(distance)[valid_samples, valid_samples]
  meta_subset <- dist_df[valid_samples, ]
  
  adonis_res <- adonis2(as.dist(dist_subset) ~ type, data = meta_subset, permutations = 999)
  print(adonis_res)
  
  cat("\n=== Pairwise Adonis ===\n")
  pairwise_res <- pairwise.adonis(as.dist(dist_subset), factors = meta_subset$type)
  print(pairwise_res)
  
  sink()
  message("  -> Saved PCoA_Statistics.txt")
  
}, error = function(e) { message("  Error in PCoA analysis: ", e$message) })

# ==============================================================================
# 4. Heatmap Visualization
# ==============================================================================
message("Generating Figure 4: Heatmap...")

tryCatch({
  heatmap_data_path <- file.path(INPUT_DIR, "core_ARG_livestock.csv")
  if (!file.exists(heatmap_data_path)) stop("Heatmap data file not found!")
  
  file <- read.csv(heatmap_data_path, row.names = 1, check.names = FALSE)
  
  pdf(file.path(OUTPUT_DIR, "Figure4_Heatmap.pdf"), width = 8, height = 10)
  pheatmap(file,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           legend_breaks = seq(-23, 10, 4),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           border_color = NA
  )
  dev.off()
  message("  -> Saved Figure4_Heatmap.pdf")
}, error = function(e) { message("  Error in Heatmap plotting: ", e$message) })

message("All visualization tasks completed successfully!")