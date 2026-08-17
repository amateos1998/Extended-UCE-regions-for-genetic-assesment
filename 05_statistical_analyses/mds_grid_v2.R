setwd("/your/path/")

library(ggplot2)
library(scales)
library(gridExtra) # For grid.arrange
library(grid) # Need this for textGrob and gpar/expression

# Function to read and prepare a single MDS plot
create_mds_plot <- function(annot_file, input_header, flip_x = FALSE, flip_y = FALSE) {
  
  # --- 1. Determine Species and Load/Filter Metadata ---
  if (grepl("atersa", input_header)) {
    species_name <- "ahya"
  } else if (grepl("aten", input_header)) {
    species_name <- "aken"
  } else if (grepl("aspat", input_header)) {
    species_name <- "aspat"
  } else {
    warning(paste("Could not determine species from header:", input_header))
    return(NULL)
  }
  
  # Read ALL metadata 
  annot_all <- read.csv(annot_file, fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  annot_species <- annot_all[annot_all$Species == species_name, ]
  
  if (nrow(annot_species) == 0) {
    warning(paste("No samples found for species", species_name, " in metadata file."))
    return(NULL)
  }
  
  # --- 2. Load MDS Coordinates ---
  mds_raw <- read.table(
    paste0(input_header, ".mds"), 
    header = TRUE, 
    stringsAsFactors = FALSE,
    row.names = 1
  )
  mds <- data.frame(Coral_ID = rownames(mds_raw), mds_raw)
  rownames(mds) <- NULL
  
  # --- 3. Outlier Removal ---
  omit_ids <- c(
    "x",
    "y"
    
  )
  
  annot_final <- annot_species[!annot_species$Coral_ID %in% omit_ids, ]
  mds_filtered <- mds[!mds$Coral_ID %in% omit_ids, ]
  
  # --- 4. ROBUST DATA ALIGNMENT ---
  PC_data <- merge(mds_filtered, annot_final, by = "Coral_ID")
  
  if (nrow(PC_data) == 0) {
    warning(paste0("Plot for ", input_header, " is EMPTY. No samples match between MDS file and filtered metadata."))
    return(NULL)
  }
  
  # --- 5. Prepare data for plotting ---
  PC <- PC_data
  colnames(PC)[2:3] <- c("PC1", "PC2") 
  
  # Extract variance explained
  tmp <- unlist(strsplit(colnames(mds_raw), split = "_"))
  val <- as.numeric(tmp[seq(2, length(tmp), 2)])
  
  # --- 6. Plotting 
  
  # Generate one unique colour per individual
  n_ind <- length(unique(PC$Coral_ID))
  ind_colours <- setNames(
    colorRampPalette(c(
      "#E41A1C","#377EB8","#4DAF4A","#984EA3",
      "#FF7F00","#A65628","#F781BF","#66C2A5",
      "#FFD700","#00CED1","#8B4513","#32CD32"
    ))(n_ind),
    unique(PC$Coral_ID)
  )
  
  p <- ggplot(PC, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = Coral_ID), size = 2, alpha = 0.7) +  # <-- was: colour = Location
    scale_color_manual(values = ind_colours) +                   # <-- was: scale_color_gradient(...)
    xlab(paste("PC1 (", ifelse(is.na(val[1]), "NA", signif(val[1], digits = 3)), "%)", sep = "")) +
    ylab(paste("PC2 (", ifelse(is.na(val[2]), "NA", signif(val[2], digits = 3)), "%)", sep = "")) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # X-axis scale 
  if (flip_x) {
    p <- p + scale_x_reverse(labels = scales::number_format(accuracy = 0.01))
  } else {
    p <- p + scale_x_continuous(labels = scales::number_format(accuracy = 0.01))
  }
  
  # Flip y-axis 
  if (flip_y) p <- p + scale_y_reverse(labels = scales::number_format(accuracy = 0.01))
  
  return(p)
}


# Function to assemble 12 MDS plots in a 3×4 grid with row and column labels
plot_mds_panel <- function(annot_file, headers, flip_x = rep(FALSE, 12), flip_y = rep(FALSE, 12), output_name = "MDS_panel") {
  if (length(headers) != 12 || length(flip_x) != 12 || length(flip_y) != 12) {
    stop("Please provide exactly 12 MDS input headers and 12-length logical vectors for flip_x and flip_y.")
  }
  
  # --- 1. Generate plots ---
  plots <- mapply(function(h, fx, fy) create_mds_plot(annot_file, h, fx, fy),
                  headers, flip_x, flip_y, SIMPLIFY = FALSE)
  
  plots <- plots[!sapply(plots, is.null)] # Remove any NULL plots
  
  if (length(plots) != 12) {
    stop("Not all 12 plots were successfully generated. Cannot proceed with fixed layout.")
  }
  
  # --- 2. Define and create text ---
  
  # A) Column Labels (4 grobs)
  col_labels_text <- c("± 500 bp", "± 5 kb", "± 20 kb", "WG")
  col_grobs <- lapply(col_labels_text, function(txt) {
    grid::textGrob(txt, gp=grid::gpar(fontsize=12, fontface="bold"))
  })
  
  # B) Row Labels (3 grobs) 
  row_labels_text <- list(
    expression(paste(italic("A. tersa"), "")), 
    expression(italic("A. kenti")), 
    expression(italic("A. spathulata"))
  )
  row_grobs <- lapply(row_labels_text, function(expr) {
    grid::textGrob(expr, gp=grid::gpar(fontsize=12, fontface="bold"), rot=90)
  })
  
  # --- 3. Combine ALL grobs into a single list (12 plots + 4 col grobs + 3 row grobs = 19 grobs) ---
  grob_list <- c(plots, col_grobs, row_grobs) 
  
  # --- 4. Define the 4x5 Grid Layout Matrix 
  
  n_plots <- 12
  n_cols_plots <- 4
  n_rows_plots <- 3 
  
  # Indices for the grob_list:
  # Plots: 1-12 
  # Col Labels: 13-16
  # Row Labels: 17-19
  
  plots_mat <- matrix(1:n_plots, nrow=n_rows_plots, ncol=n_cols_plots, byrow=TRUE)
  
  layout_mat <- matrix(NA, nrow=4, ncol=5) 
  
  layout_mat[1:3, 1] <- 17:19

  layout_mat[1:3, 2:5] <- plots_mat
 
  layout_mat[4, 2:5] <- 13:16
  

  
  # --- 5. ARRANGE AND SAVE ---
  
  
  col_widths <- c(0.08, rep(0.92/4, 4)) 
 
  row_heights <- c(rep(0.32, n_rows_plots), 0.04)
  
  gridExtra::grid.arrange(
    grobs = grob_list, 
    layout_matrix = layout_mat,
    widths = col_widths,
    heights = row_heights
  )
  
  #preview
  readline(prompt = "Press [Enter] to save the panel as SVG...")
  
  # Save as SVG 
  svg(filename = paste0(output_name, ".svg"), width = 9, height = 7) 
  gridExtra::grid.arrange(
    grobs = grob_list, 
    layout_matrix = layout_mat,
    widths = col_widths,
    heights = row_heights
  )
  dev.off()
  
  message("Panel saved as ", output_name, ".svg")
}

# ===============================
#USE
# ===============================

headers <- c(
  
  "top100_60Samples_to_atersa_500_MinQ20_SNPpval1e-3_Prior1", "top100_60Samples_to_atersa_5k_MinQ20_SNPpval1e-3_Prior1",
  "top100_60Samples_to_atersa_20k_MinQ20_SNPpval1e-3_Prior1", "60Samples_to_atersa_WG_MinQ20_SNPpval1e-3_Prior1",
  
  "top100_71Samples_to_aten_500_MinQ20_SNPpval1e-3_Prior1", "top100_71Samples_to_aten_5k_MinQ20_SNPpval1e-3_Prior1",
  "top100_71Samples_to_aten_20k_MinQ20_SNPpval1e-3_Prior1","71Samples_to_aten_WG_MinQ20_SNPpval1e-3_Prior1",
  
  "top100_78Samples_to_aspat_500_MinQ20_SNPpval1e-3_Prior1", "top100_78Samples_to_aspat_5k_MinQ20_SNPpval1e-3_Prior1", 
  "top100_78Samples_to_aspat_20k_MinQ20_SNPpval1e-3_Prior1","78Samples_to_aspat_WG_MinQ20_SNPpval1e-3_Prior1"
)

# Flip X-axis
flip_x <- rep(FALSE, 12)
flip_x[c(6,7)] <- TRUE

# Flip y-axis 
flip_y <- rep(FALSE, 12)
flip_y[c(5,7)] <- TRUE

# Use
plot_mds_panel("metadata_all.txt", headers, flip_x, flip_y, output_name = "SUPP_mainpops_mds_top100_loci")


