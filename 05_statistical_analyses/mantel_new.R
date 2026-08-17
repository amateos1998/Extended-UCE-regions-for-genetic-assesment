library(vegan)
library(ggplot2)
library(reshape2)
library(patchwork)
library(ggrepel)

prefix <- "aten_main_top200"
setwd("C:/Users/amateos/OneDrive - Australian Institute of Marine Science/AIMS_Semesters 3 and 4/Semesters 3 and 4/backups_from_HPC/results_Mantel")
# 1. INITIAL MANTEL TEST (no filtering, to id outliers)
{
  file_names <- paste0(prefix, c("_500", "_5k","_20k", "_WG"))
  
  
  matrices <- list()
  
  for (name in file_names) {
    file_path <- paste0(name, ".dist")
    mat <- read.table(file_path, sep = "\t", skip = 2, row.names = 1, check.names = FALSE)
    matrices[[name]] <- as.dist(as.matrix(mat))
  }
  
  names(matrices) <- file_names
  n_files <- length(file_names)
  results <- matrix(NA, nrow = length(file_names), ncol = length(file_names), dimnames = list(file_names, file_names))
  mantel_data <- data.frame()
  
  for (i in 1:n_files) {
    for (j in i:n_files) {
      if (i != j) {
        test <- mantel(matrices[[i]], matrices[[j]], method = "pearson", permutations = 9999)
        results[i, j] <- test$statistic
        results[j, i] <- test$signif
        
        mantel_data <- rbind(mantel_data, data.frame(
          Matrix1 = file_names[i],
          Matrix2 = file_names[j],
          r = round(test$statistic, 3),
          p = ifelse(test$signif == 0, "≤ 0.001", as.character(round(test$signif, 3)))
        ))
      }
    }
  }
}

# 2. INDIVIDUAL SCATTERPLOTS vs. WG + IDENTIFY OUTLIERS
{
  wg_matrix <- as.matrix(as.dist(matrices[[paste0(prefix, "_WG")]]))
  outlier_pairs <- c()
  
  for (name in file_names) {
    if (name != paste0(prefix, "_WG")) {
      mat_i <- as.matrix(as.dist(matrices[[name]]))
      common_samples <- intersect(rownames(mat_i), rownames(wg_matrix))
      mat_i <- mat_i[common_samples, common_samples]
      mat_wg <- wg_matrix[common_samples, common_samples]
      
      pair_indices <- which(lower.tri(mat_i), arr.ind = TRUE)
      labels <- paste(rownames(mat_i)[pair_indices[,1]], rownames(mat_i)[pair_indices[,2]], sep = " vs ")
      
      df <- data.frame(
        x = mat_i[lower.tri(mat_i)],
        y = mat_wg[lower.tri(mat_wg)],
        pair = labels
      )
      
      # Identify outliers
      df$label <- ifelse(df$x < 0.14 | df$y < 0.100, df$pair, "")
      outlier_pairs <- c(outlier_pairs, df$label[df$label != ""])
      
      test <- mantel(as.dist(mat_i), as.dist(mat_wg), method = "pearson", permutations = 999)
      r_val <- round(test$statistic, 3)
      p_val <- ifelse(test$signif == 0, "≤ 0.001", as.character(round(test$signif, 3)))
      
      p <- ggplot(df, aes(x = x, y = y)) +
        geom_point(size = 0.8, alpha = 0.9, colour = "black") +
        ggrepel::geom_text_repel(aes(label = label), size = 2.5, max.overlaps = Inf,
                                 box.padding = 0.5, point.padding = 0.5, force = 2, force_pull = 0.1) +
        geom_smooth(method = "lm", se = FALSE, colour = "red") +
        annotate("text", x = min(df$x), y = max(df$y), hjust = 0, vjust = 1,
                 label = paste0("r = ", r_val, "\np = ", p_val), size = 4) +
        labs(x = NULL, y = NULL, title = NULL) +
        theme_minimal() +
        theme(
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank()
        )
      
      print(p)
    }
  }
  
  # Extract unique sample names from outlier pairs
  outlier_samples <- unique(unlist(strsplit(outlier_pairs, " vs ")))
  cat("\nOutlier samples identified (", length(outlier_samples), "):\n",
      paste(outlier_samples, collapse = ", "), "\n")
}



##########EVERYTHING BELOW IS OPTIONAL##############


# 3. RE-RUN MANTEL TEST AFTER REMOVING OUTLIERS (only if outliers exist)
if (length(outlier_samples) > 0) {
  
  matrices <- list()
  
  for (name in file_names) {
    file_path <- paste0(name, ".dist")
    mat <- read.table(file_path, sep = "\t", skip = 2, row.names = 1, check.names = FALSE)
    
    # Fix shape and align names
    mat <- as.matrix(mat)
    if (ncol(mat) > nrow(mat)) {
      mat <- mat[, 1:nrow(mat)]
    }
    colnames(mat) <- rownames(mat)
    
    shared <- setdiff(rownames(mat), outlier_samples)
    
    if (length(shared) >= 2) {
      mat <- mat[shared, shared]
      matrices[[name]] <- as.dist(mat)
      cat(name, ": kept", length(shared), "samples\n")
    } else {
      warning(paste(name, "- skipped: too few samples after filtering"))
    }
  }
  
  results <- matrix(NA, nrow = 7, ncol = 7, dimnames = list(file_names, file_names))
  mantel_data <- data.frame()
  
  for (i in 1:7) {
    for (j in i:7) {
      if (i != j && !is.null(matrices[[i]]) && !is.null(matrices[[j]])) {
        test <- mantel(matrices[[i]], matrices[[j]], method = "pearson", permutations = 999)
        results[i, j] <- test$statistic
        results[j, i] <- test$signif
      }
    }
  }
  
  # 3.5 SCATTERPLOTS vs. WG AFTER REMOVAL
  wg_matrix <- as.matrix(as.dist(matrices[["aspat_WG"]]))
  
  for (name in file_names) {
    if (name != paste0(prefix, "_WG") && !is.null(matrices[[name]])) {
      mat_i <- as.matrix(as.dist(matrices[[name]]))
      
      common_samples <- intersect(rownames(mat_i), rownames(wg_matrix))
      mat_i <- mat_i[common_samples, common_samples]
      mat_wg <- wg_matrix[common_samples, common_samples]
      
      df <- data.frame(
        x = mat_i[lower.tri(mat_i)],
        y = mat_wg[lower.tri(mat_wg)]
      )
      
      test <- mantel(as.dist(mat_i), as.dist(mat_wg), method = "pearson", permutations = 999)
      r_val <- round(test$statistic, 3)
      p_val <- ifelse(test$signif == 0, "≤ 0.001", as.character(round(test$signif, 3)))
      
      p <- ggplot(df, aes(x = x, y = y)) +
        geom_point(size = 0.8, alpha = 0.9, colour = "blue") +
        geom_smooth(method = "lm", se = FALSE, colour = "red") +
        annotate("text", x = min(df$x), y = max(df$y), hjust = 0, vjust = 1,
                 label = paste0("r = ", r_val, "\np = ", p_val), size = 4) +
        labs(title = paste(name, "vs.", paste0(prefix, "_WG"), "(after filtering)")) +
        theme_minimal() +
        theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
      
      print(p)
    }
  }
  
} else {
  cat("No outliers detected — skipping Steps 3 and 3.5.\n")
}



# 4. HEATMAP OF FILTERED MANTEL RESULTS (after outlier removal)
{
  heatmap_df <- melt(results, varnames = c("Matrix1", "Matrix2"), value.name = "Value")
  
  gg_heatmap <- ggplot(heatmap_df, aes(x = Matrix1, y = Matrix2, fill = Value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.3f", Value)), size = 3) +
    scale_fill_gradient2(
      low = "steelblue3",
      mid = "#FFF7BC",
      high = "red",
      midpoint = 0.5,
      name = "r / p values",
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1"),
      limits = c(0, 1),
      na.value = "lightgrey"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      legend.key.height = unit(1.5, "cm")
    )
  
  print(gg_heatmap)
  ggsave(paste0(prefix, "_mantel_heatmap_filtered.png"), gg_heatmap, width = 8, height = 6, dpi = 300)
}


# 5. SCATTERPLOT MATRIX OF BASE MATRICES vs WG
{
  base_names <- paste0(prefix, c("_500", "_5k", "_20k"))
  wg_name <- paste0(prefix, "_WG")
  

  
  plot_list <- lapply(base_names, function(base) {
    x <- as.vector(matrices[[base]])
    y <- as.vector(matrices[[wg_name]])
    df <- data.frame(x = x, y = y)
    
    test <- mantel(matrices[[base]], matrices[[wg_name]], method = "pearson", permutations = 999)
    r_val <- round(test$statistic, 3)
    p_val <- ifelse(test$signif == 0, "≤ 0.001", as.character(round(test$signif, 3)))
    
    # Get plot-specific min/max for annotation placement (since global limits are gone)
    x_min <- min(df$x, na.rm = TRUE)
    y_max <- max(df$y, na.rm = TRUE)
    
    ggplot(df, aes(x = x, y = y)) +
      geom_point(size = 0.8, alpha = 0.9, colour = "black") +
      geom_smooth(method = "lm", se = FALSE, colour = "red") +
      # Annotate position adjusted back to use plot-specific limits
      annotate("text", x = x_min, y = y_max, hjust = 0, vjust = 1, 
               label = paste0("r = ", r_val, "\np = ", p_val), size = 4) +
      theme_minimal(base_size = 13) +
      
      labs(title = NULL, x = NULL, y = NULL) +
      # *** MODIFIED: Keep coord_fixed() for 1:1 aspect ratio, but remove xlim/ylim ***
      coord_fixed() +
      theme(
        # Zero margin for each individual plot
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
      )
  })
  
  # Arrange all plots in a single horizontal row
  scatter_matrix <- wrap_plots(plot_list, nrow = 1, widths = c(1, 1, 1), heights = 1) 
    # Apply zero spacing/margins to the final patchwork layout
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          panel.spacing = unit(0, "cm"))
  
  print(scatter_matrix)
  
  ggsave(paste0(prefix, "_mantel_scatter_matrix_horizontal_no_labels_grid.svg"),
         scatter_matrix, width = 7.5, height = 2.5, units = "cm")
  #ggsave(paste0(prefix, "_mantel_scatter_matrix_horizontal_no_labels_grid.svg"),
         #scatter_matrix, width = 12, height = 6)
}

