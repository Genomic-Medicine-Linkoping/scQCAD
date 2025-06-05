# Helper function to check for the presence of clonotype columns
check_clonotype <- function(seurat_obj, cell_type) {
  clonotype_col <- paste0(cell_type, "_clonotype_id")
  if (clonotype_col %in% colnames(seurat_obj@meta.data)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Clonotype Frequency Analysis
clonotype_freq_analysis <- function(seurat_obj, cell_type, seurat_out_dir) {
  clonotype_col <- paste0(cell_type, "_clonotype_id")
  cat("Performing clonotype frequency analysis for", cell_type, "cells...\n")

  if (check_clonotype(seurat_obj, cell_type)) {
    # Frequency table
    clonotype_freq <- table(seurat_obj@meta.data[[clonotype_col]])
    print(clonotype_freq)

    # Plot clonotype frequency
    clonotype_freq_df <- as.data.frame(clonotype_freq)
    freq_histogram_vdj <- ggplot2::ggplot(
      clonotype_freq_df, aes(x = reorder(Var1, -Freq), y = Freq)
    ) +
      ggplot2::geom_bar(stat = "identity") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1)
      ) +
      xlab("Clonotype ID") +
      ylab("Frequency") +
      ggtitle(paste(cell_type, "Clonotype Frequency"))

    filename <- paste0(
      seurat_out_dir, "/frequency_vdj_", cell_type, "_cell.jpeg"
    )
    ggplot2::ggsave(
      filename = filename,
      plot = freq_histogram_vdj
    )
  }
}


# UMAP Visualization of Clonotypes
clonotype_umap_plot <- function(seurat_obj, cell_type, seurat_out_dir) {
  clonotype_col <- paste0(cell_type, "_clonotype_id")
  cat("Generating UMAP plot colored by", cell_type, "clonotypes...\n")
  if (check_clonotype(seurat_obj, cell_type)) {
    clonotype_umap <- DimPlot(
      seurat_obj,
      group.by = clonotype_col,
      reduction = "umap",
      label = TRUE,
      repel = TRUE
    ) + NoLegend() +
      ggtitle(paste("UMAP Colored by", cell_type, "Clonotype"))

    filename <- paste0(seurat_out_dir, "/umap_vdj_", cell_type, "_cell.jpeg")
    ggplot2::ggsave(
      filename = filename,
      plot = clonotype_umap
    )
  }
}

# Highlight cells with a specific clonotype on UMAP
highlight_clonotype_umap <- function(
    seurat_obj, cell_type, clonotype_id, seurat_out_dir) {
  clonotype_col <- paste0(cell_type, "_clonotype_id")
  cat(
    "Highlighting cells with", cell_type, "clonotype",
    clonotype_id, "on UMAP...\n"
  )
  if (check_clonotype(seurat_obj, cell_type)) {
    highlighted_clonotype_umap <- DimPlot(
      seurat_obj,
      cells.highlight = WhichCells(
        seurat_obj,
        expression = get(clonotype_col) == clonotype_id
      )
    ) +
      ggtitle(paste("Cells with", cell_type, "Clonotype", clonotype_id))
    filename <- paste0(
      seurat_out_dir, "/highlighted_umap_vdj_", cell_type, "_cell.jpeg"
    )
    ggplot2::ggsave(
      filename = filename,
      plot = highlighted_clonotype_umap
    )
  }
}

# Differential Expression Analysis between clonotype groups
clonotype_diff_expr <- function(
    seurat_obj, cell_type, seurat_out_dir) {
  clonotype_col <- paste0(cell_type, "_clonotype_id")
  cat(
    "Performing differential expression analysis for",
    cell_type, "clonotype ...\n"
  )
  if (check_clonotype(seurat_obj, cell_type)) {
    for (curr_clonotype_id in unique(seurat_obj@meta.data[[clonotype_col]])) {
      if (!is.na(curr_clonotype_id)) {
        print(clonotype_col)
        # Get the cells corresponding to the current clonotype
        clonotype_cells <- seurat_obj@meta.data[[clonotype_col]] == curr_clonotype_id
        cn <- rownames(seurat_obj@meta.data)[which(clonotype_cells)]
        # Check if the current clonotype has at least 3 cells
        if (length(cn) < 3) {
          cat(
            "Skipping clonotype", curr_clonotype_id,
            "because it has fewer than 3 cells.\n"
          )
          next
        }
        print(curr_clonotype_id)

        de_results <- FindMarkers(
          seurat_obj,
          ident.1 = curr_clonotype_id,
          group.by = clonotype_col
        )
        de_results[["Gene"]] <- row.names(de_results)
        ordered_cols <- c(
          "Gene",
          names(de_results)[names(de_results) != "Gene"]
        )
        de_results <- de_results %>%
          dplyr::select(all_of(ordered_cols)) %>%
          dplyr::mutate(
            across(
              where(is.numeric),
              ~ ifelse(. < 0.001,
                formatC(., format = "e", digits = 2), round(., 3)
              )
            )
          )


        print(head(de_results))
        write.csv(
          de_results,
          file = paste0(
            seurat_out_dir, "/diff_expr_vdj_", cell_type,
            "_cell_", curr_clonotype_id, "_vs_all.csv"
          ),
          row.names = FALSE
        )
      }
    }
  }
}


vdj_downstream_analysis <- function(
    seurat_obj, cell_type = c(), seurat_out_dir) {
  # Perform analyses based on the input cell type
  if ("t" %in% cell_type) {
    # T-cell VDJ analyses
    clonotype_freq_analysis(seurat_obj, "t", seurat_out_dir)
    clonotype_umap_plot(seurat_obj, "t", seurat_out_dir)
    clonotype_diff_expr(seurat_obj, "t", seurat_out_dir)
  }

  if ("b" %in% cell_type) {
    # B-cell VDJ analyses
    clonotype_freq_analysis(seurat_obj, "b", seurat_out_dir)
    clonotype_umap_plot(seurat_obj, "b", seurat_out_dir)
    clonotype_diff_expr(seurat_obj, "b", seurat_out_dir)
  }
}
