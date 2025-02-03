load_data_10X <- function(data_dir = NULL, data_file = NULL) {
  if (!is.null(data_dir)) {
    if (dir.exists(data_dir)) {
      data_10_x <- Seurat::Read10X(data.dir = data_dir)
      return(data_10_x)
    }
  } else if (!is.null(data_file)) {
    if (file.exists(data_file)) {
      data_10_x <- Seurat::Read10X_h5(filename = data_file)
      return(data_10_x)
    }
  } else {
    cat("10X data not found!")
  }
}

# Function to load metadata and integrate cell data
metadata <- function(data_10_x, aggregate_csv_file) {
  if (!is.null(aggregate_csv_file)) {
    if (file.exists(aggregate_csv_file)) {
      samples <- read.csv(aggregate_csv_file, stringsAsFactors = FALSE)

      if (is.list(data_10_x)) {
        raw_data_10_x <- data_10_x[["Gene Expression"]]
      } else {
        raw_data_10_x <- data_10_x
      }

      cells <- new("seurat", raw.data = raw_data_10_x)

      cellcodes <- as.data.frame(cells@raw.data@Dimnames[[2]])
      colnames(cellcodes) <- "barcodes"
      rownames(cellcodes) <- cellcodes$barcodes
      cellcodes$libcodes <- as.factor(
        gsub(pattern = ".+-", replacement = "", cellcodes$barcodes)
      )
      for (mcol in names(samples)) {
        if (mcol != "sample_outs") {
          cellcodes[[mcol]] <- as.vector(samples[[mcol]][cellcodes$libcodes])
        }
      }

      return(cellcodes)
    }
  } else {
    return(NULL)
  }
}

seurat_set_orig_ident <- function(seurat_obj, sample_identity) {
  meta_cols <- names(sample_identity)
  if ("sample_id" %in% meta_cols) {
    seurat_obj@meta.data["orig.ident"] <- seurat_obj@meta.data["sample_id"]
    return(seurat_obj)
  } else {
    seurat_obj@meta.data["orig.ident"] <- gsub(".+-", "", rownames(seurat_obj))
    return(seurat_obj)
  }
}

acronym <- function(name_list) {
  acronym_list <- lapply(
    lapply(
      sapply(name_list, function(x) stringr::str_split(x, " ")),
      function(y) stringr::str_sub(y, 1, 1)
    ),
    function(w) paste0(w, collapse = "")
  )
  return(acronym_list)
}
# Generic function to create a Seurat object from Gene Expression and
# add additional assays
create_seurat_object <- function(
    data_10_x, sample_identity, project_name) {
  # Create base Seurat object using Gene Expression

  if (is.list(data_10_x)) {
    seurat_obj <- Seurat::CreateSeuratObject(
      counts = data_10_x[["Gene Expression"]],
      meta.data = sample_identity,
      project = project_name
    )

    # Check if additional assays are present and add them
    additional_assays <- setdiff(names(data_10_x), "Gene Expression")
    acronym_list <- acronym(additional_assays)
    for (assay in additional_assays) {
      current_assay <- SeuratObject::CreateAssay5Object(
        counts = data_10_x[[assay]]
      )
      seurat_obj[[acronym_list[assay][[1]]]] <- current_assay
    }
  } else {
    seurat_obj <- Seurat::CreateSeuratObject(
      counts = data_10_x,
      meta.data = sample_identity,
      project = project_name
    )
  }

  return(seurat_obj)
}

add_clonotype <- function(tcr_file, seurat_obj, type = "t") {
  if (is.null(tcr_file)) {
    return(seurat_obj)
  }

  tcr_prefix <- dirname(normalizePath(tcr_file))
  tcr <- read.csv(
    paste(tcr_prefix, "filtered_contig_annotations.csv", sep = "/")
  )

  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  tcr <- tcr[!duplicated(tcr$barcode), ]

  # Only keep the barcode and clonotype columns.
  # We'll get additional clonotype info from the clonotype table.
  tcr <- tcr[, c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_prefix, "clonotypes.csv", sep = "/"))

  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"

  # Reorder so barcodes are first column and set them as rownames.
  tcr <- tcr[, c(2, 1, 3)]
  rownames(tcr) <- tcr[, 1]
  tcr[, 1] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep = "_")
  # Add to the Seurat object's metadata.
  seurat_obj <- SeuratObject::AddMetaData(object = seurat_obj, metadata = tcr)
  return(seurat_obj)
}

seurat_qc <- function(seurat_obj) {
  # Add mitochondrial and ribosomal percentages
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-"
  )
  seurat_obj[["percent.rb"]] <- Seurat::PercentageFeatureSet(
    seurat_obj,
    pattern = "^RP[SL]"
  )
  return(seurat_obj)
}

quantiles_for_qc <- function(seurat_obj) {
  probs <- c(0.90, 0.95, 0.99)
  # Calculate quantiles for all features of interest
  nFeature_quantiles <- quantile(seurat_obj$nFeature_RNA, probs = probs)
  nCount_quantiles <- quantile(seurat_obj$nCount_RNA, probs = probs)
  mt_quantiles <- quantile(seurat_obj$percent.mt, probs = probs)
  rb_quantiles <- quantile(seurat_obj$percent.rb, , probs = probs)

  # Round quantiles for better readability
  nFeature_quantiles <- round(nFeature_quantiles, 0)
  nCount_quantiles <- round(nCount_quantiles, 0)
  mt_quantiles <- round(mt_quantiles, 2)
  rb_quantiles <- round(rb_quantiles, 2)

  return(
    list(
      nFeature_quantiles,
      nCount_quantiles,
      mt_quantiles,
      rb_quantiles
    )
  )
}

seurat_filtering <- function(
    seurat_obj,
    min_cells,
    min_features,
    max_features,
    percent_mt,
    percent_rb) {
  quantile_list <- quantiles_for_qc(seurat_obj)

  if (is.null(percent_mt)) {
    percent_mt <- quantile_list[[3]][2]
  }
  if (is.null(percent_rb)) {
    percent_rb <- quantile_list[[4]][2]
  }

  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > min_features & nFeature_RNA < max_features &
      nCount_RNA > min_cells & percent.mt < percent_mt & percent.rb < percent_rb
  )
  return(seurat_obj)
}

split_seurat_by_layer <- function(
    seurat_obj, assay = "RNA", layer_column = NULL) {
  if (!is.null(layer_column)) {
    seurat_obj[["RNA"]] <- S4Vectors::split(
      seurat_obj[["RNA"]],
      f = seurat_obj@meta.data[[layer_column]]
    )
  }
  return(seurat_obj)
}

method_dict <- list(
  "CCAIntegration" = c("pca", "integrated.cca"),
  "RPCAIntegration" = c("pca", "integrated.rpca"),
  "HarmonyIntegration" = c("pca", "harmony"),
  "JointPCAIntegration" = c("pca", "integrated.dr")
)

# Function to integrate Seurat layers using specified method
integrate_seurat_layers <- function(
    seurat_obj, integration_method = NULL, enable_sct = TRUE) {
  if (is.null(integration_method)) {
    return(seurat_obj)
  } else {
    if (enable_sct) {
      normalization_method <- "SCT"
    } else {
      normalization_method <- "LogNormalize"
    }

    seurat_obj <- Seurat::IntegrateLayers(
      seurat_obj,
      method = integration_method,
      orig.reduction = method_dict[[integration_method]][1],
      new.reduction = method_dict[[integration_method]][2],
      normalization.method = normalization_method
    )
    return(seurat_obj)
  }
}

join_seurat_layers <- function(seurat_obj, assay = "RNA", layer_column = NULL) {
  if (!is.null(layer_column)) {
    seurat_obj[[assay]] <- SeuratObject::JoinLayers(
      seurat_obj[[assay]]
    )
  }
  return(seurat_obj)
}



seurat_normalization <- function(seurat_obj, sct = FALSE) {
  if (sct) {
    seurat_obj <- Seurat::SCTransform(
      seurat_obj,
      vars.to.regress = "percent.mt",
      verbose = TRUE
    )

    seurat_obj <- Seurat::PrepSCTFindMarkers(seurat_obj)

    return(seurat_obj)
  } else {
    seurat_obj <- Seurat::NormalizeData(
      seurat_obj,
      normalization.method = "LogNormalize",
      scale.factor = 10000
    )
    seurat_obj <- Seurat::FindVariableFeatures(
      seurat_obj,
      selection.method = "vst", nfeatures = 2000
    )
    seurat_obj <- Seurat::ScaleData(
      seurat_obj,
      features = rownames(seurat_obj),
      vars.to.regress = "percent.mt"
    )
    return(seurat_obj)
  }
}

seurat_dim_reduction_linear <- function(seurat_obj) {
  n_cells <- ncol(seurat_obj)
  npcs <- ifelse(n_cells < 50, min(n_cells, nrow(seurat_obj)) - 1, 50)

  seurat_obj <- Seurat::RunPCA(
    seurat_obj,
    npcs = npcs,
    verbose = TRUE
  )
  return(seurat_obj)
}

seurat_dim_reduction_nonlinear <- function(
    seurat_obj, integration_method = NULL) {
  if (is.null(integration_method)) {
    reduction <- "pca"
    reduction_name <- "umap"
  } else {
    reduction <- method_dict[[integration_method]][2]
    reduction_name <- paste0("umap.", reduction)
  }
  seurat_obj <- Seurat::RunUMAP(
    seurat_obj,
    reduction = reduction,
    dims = 1:30,
    verbose = TRUE,
    reduction.name = reduction_name
  )
  return(seurat_obj)
}


seurat_clustering <- function(
    seurat_obj, integration_method = NULL, resolution = 0.5) {
  if (is.null(integration_method)) {
    reduction <- "pca"
  } else {
    reduction <- method_dict[[integration_method]][2]
  }
  seurat_obj <- Seurat::FindNeighbors(
    seurat_obj,
    reduction = reduction,
    dims = 1:30
  )
  seurat_obj <- Seurat::FindClusters(
    seurat_obj,
    resolution = resolution,
    random.seed = 1234
  )
}

seurat_annotation <- function(
    seurat_obj, species = "human", ref.data = NULL, ensembl = FALSE, num_threads = 8) {
  cat("Using singleR for annotation ...")

  if (is.null(ref.data)) {
    if (species == "human") {
      ref.data <- celldex::HumanPrimaryCellAtlasData(ensembl)
    }
    if (species == "mouse") {
      ref.data <- celldex::MouseRNAseqData(ensembl)
    }
  }

  if (is.null(ref.data)) {
    stop(paste0("Error: ", species, " does not have a reference available."))
  }

  sc_obj <- Seurat::as.SingleCellExperiment(seurat_obj)
  # Performing predictions.
  predictions <- SingleR::SingleR(
    test = sc_obj,
    ref = ref.data,
    labels = ref.data$label.main,
    num.threads = num_threads
  )

  print(table(predictions$labels))

  seurat_obj[["SingleR.labels"]] <- predictions$labels

  return(seurat_obj)
}

seurat_all_markers <- function(
    seurat_obj,
    assay = "RNA",
    cluster_column = "seurat_clusters",
    prefix, seurat_out_dir) {
  Idents(object = seurat_obj) <- seurat_obj@meta.data[[cluster_column]]

  all_markers <- Seurat::FindAllMarkers(
    seurat_obj,
    assay = assay,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

  write.csv(
    all_markers,
    file = paste0(
      seurat_out_dir,
      "/",
      prefix,
      "_",
      cluster_column,
      "_all_markers.csv"
    ),
    row.names = FALSE
  )

  seurat_top <- function(df, n) {
    return(
      df %>%
        dplyr::group_by(cluster) %>%
        dplyr::arrange(avg_log2FC, .by_group = TRUE) %>%
        dplyr::slice_head(n = n) %>%
        dplyr::ungroup()
    )
  }
  top1 <- seurat_top(all_markers, 1)
  top6 <- seurat_top(all_markers, 6)
  top10 <- seurat_top(all_markers, 10)
  l <- list(all_markers, top1, top6, top10)
  return(l)
}

# Perform DE analysis within the same cell type across conditions
seurat_de <- function(
    seurat_obj,
    annotation_column = NULL,
    assay = "RNA",
    seurat_out_dir) {
  original_idents <- get_idents_column_name(seurat_obj)

  SeuratObject::Idents(seurat_obj) <- annotation_column

  celltypes <- unique(seurat_obj@meta.data[[annotation_column]])
  celltypepairs <- combn(celltypes, 2)

  for (i in 1:ncol(celltypepairs)) {
    pair <- celltypepairs[, i]

    ident_1 <- pair[1]
    ident_2 <- pair[2]

    # Get the cells corresponding to ident_1 and ident_2
    cells_1 <- SeuratObject::WhichCells(seurat_obj, idents = ident_1)
    cells_2 <- SeuratObject::WhichCells(seurat_obj, idents = ident_2)

    # Check if both groups have at least 3 cells
    if (!(length(cells_1) >= 3 & length(cells_2) >= 3)) {
      next
    }


    curr_celltype_de <- Seurat::FindMarkers(
      seurat_obj,
      assay = assay,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      ident.1 = ident_1,
      ident.2 = ident_2
    )

    curr_celltype_de["Comparison"] <- paste0(ident_1, "_Vs_", ident_2)
    curr_celltype_de["Gene"] <- rownames(curr_celltype_de)

    rownames(curr_celltype_de) <- c(1:nrow(curr_celltype_de))

    curr_celltype_de <- curr_celltype_de %>%
      dplyr::select(
        all_of(
          c("Gene", "Comparison", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
        )
      )

    if (i == 1) {
      all_celltype_de <- curr_celltype_de
    } else {
      all_celltype_de <- rbind(all_celltype_de, curr_celltype_de)
    }
  }

  SeuratObject::Idents(seurat_obj) <- original_idents

  if (exists("all_celltype_de")) {
    file2save <- paste0(
      seurat_out_dir, "/de_for_", annotation_column, "_annotations.csv"
    )

    write.csv(
      all_celltype_de,
      file = file2save,
      row.names = FALSE
    )
  }

  plot_de_genes(
    seurat_obj,
    all_celltype_de,
    annotation_column,
    NULL,
    seurat_out_dir,
    paste0("de_between_", annotation_column, "_levels"),
    assay
  )
}


# Perform DE analysis within the same cell type across conditions
seurat_de_across_condition_same_cell_type <- function(
    seurat_obj,
    annotation_column = NULL,
    condition_column = NULL,
    assay = "RNA",
    seurat_out_dir) {
  original_idents <- get_idents_column_name(seurat_obj)

  annotation_list <- seurat_obj[[annotation_column]][[1]]
  condition_list <- seurat_obj[[condition_column]][[1]]

  seurat_obj$celltype.condition <- comprehenr::to_vec(
    for (i in seq_along(annotation_list)) {
      paste0(annotation_list[i], "_", condition_list[i])
    }
  )

  SeuratObject::Idents(seurat_obj) <- "celltype.condition"

  celltypes <- unique(seurat_obj@meta.data[[annotation_column]])
  conditions <- levels(factor(unique(seurat_obj@meta.data[[condition_column]])))

  ident_levels <- levels(SeuratObject::Idents(seurat_obj))

  for (icelltype in seq_along(celltypes)) {
    celltype <- celltypes[icelltype]

    ident_1 <- paste0(celltype, "_", conditions[1])
    ident_2 <- paste0(celltype, "_", conditions[2])

    if ((!(ident_1 %in% ident_levels)) || (!(ident_2 %in% ident_levels))) {
      next
    }

    # Get the cells corresponding to ident_1 and ident_2
    cells_1 <- SeuratObject::WhichCells(seurat_obj, idents = ident_1)
    cells_2 <- SeuratObject::WhichCells(seurat_obj, idents = ident_2)

    # Check if both groups have at least 3 cells
    if (!(length(cells_1) >= 3 & length(cells_2) >= 3)) {
      next
    }


    curr_celltype_de <- Seurat::FindMarkers(
      seurat_obj,
      assay = assay,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      ident.1 = ident_1,
      ident.2 = ident_2
    )

    curr_celltype_de["Celltype"] <- celltype
    curr_celltype_de["Comparison"] <- paste0(conditions[1], "_Vs_", conditions[2])
    curr_celltype_de["Gene"] <- rownames(curr_celltype_de)

    rownames(curr_celltype_de) <- c(1:nrow(curr_celltype_de))

    curr_celltype_de <- curr_celltype_de %>%
      dplyr::select(
        all_of(
          c("Gene", "Celltype", "Comparison", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
        )
      )

    if (icelltype == 1) {
      all_celltype_de <- curr_celltype_de
    } else {
      all_celltype_de <- rbind(all_celltype_de, curr_celltype_de)
    }
  }

  SeuratObject::Idents(seurat_obj) <- original_idents

  if (exists("all_celltype_de")) {
    file2save <- paste0(
      seurat_out_dir, "/de_same_cell_for_", annotation_column, "_across_", condition_column, ".csv"
    )

    write.csv(
      all_celltype_de,
      file = file2save,
      row.names = FALSE
    )

    plot_de_genes(
      seurat_obj,
      all_celltype_de,
      annotation_column,
      condition_column,
      seurat_out_dir,
      paste0("de_within_", annotation_column, "_levels_across_", condition_column),
      assay
    )
  }
}


# Perform DE analysis with pseudo bulking
seurat_de_pseudo_bulk_across_condition <- function(
    seurat_obj,
    annotation_column = NULL,
    condition_column = NULL,
    group_by = NULL,
    assay = "RNA",
    seurat_out_dir) {
  # pseudobulk the counts based on donor-condition-celltype
  pseudo_seurat_obj <- AggregateExpression(
    seurat_obj,
    assays = assay,
    return.seurat = T,
    group.by = group_by
  )

  annotation_list <- pseudo_seurat_obj[[annotation_column]][[1]]
  condition_list <- pseudo_seurat_obj[[condition_column]][[1]]

  pseudo_seurat_obj$celltype.condition <- comprehenr::to_vec(
    for (i in seq_along(condition_list)) {
      paste0(condition_list[i], "_", annotation_list[i])
    }
  )

  SeuratObject::Idents(pseudo_seurat_obj) <- "celltype.condition"

  ident_levels <- levels(SeuratObject::Idents(pseudo_seurat_obj))

  celltypes <- unique(pseudo_seurat_obj@meta.data[[annotation_column]])
  conditions <- levels(factor(unique(pseudo_seurat_obj@meta.data[[condition_column]])))

  for (i in seq_along(celltypes)) {
    celltype <- celltypes[i]

    ident_1 <- paste0(conditions[1], "_", celltype)
    ident_2 <- paste0(conditions[2], "_", celltype)

    if ((!(ident_1 %in% ident_levels)) || (!(ident_2 %in% ident_levels))) {
      next
    }

    # Get the cells corresponding to ident_1 and ident_2
    cells_1 <- SeuratObject::WhichCells(
      pseudo_seurat_obj,
      idents = ident_1
    )
    cells_2 <- SeuratObject::WhichCells(
      pseudo_seurat_obj,
      idents = ident_2
    )

    # Check if both groups have at least 3 cells
    if (!(length(cells_1) >= 3 & length(cells_2) >= 3)) {
      next
    }

    curr_celltype_de_pseudobulk <- Seurat::FindMarkers(
      pseudo_seurat_obj,
      assay = assay,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      ident.1 = ident_1,
      ident.2 = ident_2,
      # min.cells.group = 1,
      test.use = "DESeq2"
    )

    curr_celltype_de_pseudobulk["Celltype"] <- celltype
    curr_celltype_de_pseudobulk["Comparison"] <- paste0(conditions[1], "_Vs_", conditions[2])
    curr_celltype_de_pseudobulk["Gene"] <- rownames(curr_celltype_de_pseudobulk)

    rownames(curr_celltype_de_pseudobulk) <- c(1:nrow(curr_celltype_de_pseudobulk))

    curr_celltype_de_pseudobulk <- curr_celltype_de_pseudobulk %>%
      dplyr::select(
        all_of(
          c("Gene", "Celltype", "Comparison", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
        )
      )

    if (i == 1) {
      all_celltype_de_pseudo_bulk <- curr_celltype_de_pseudobulk
    } else {
      all_celltype_de_pseudo_bulk <- rbind(all_celltype_de_pseudo_bulk, curr_celltype_de_pseudobulk)
    }
  }


  if (exists("all_celltype_de_pseudo_bulk")) {
    file2save <- paste0(
      seurat_out_dir, "/de_pseudo_bulk_for_", annotation_column, "_across_", condition_column, ".csv"
    )

    write.csv(
      all_celltype_de_pseudo_bulk,
      file = file2save,
      row.names = FALSE
    )

    plot_de_genes(
      seurat_obj,
      all_celltype_de_pseudo_bulk,
      annotation_column,
      condition_column,
      seurat_out_dir,
      paste0("pseudo_bulk_de_within_", annotation_column, "_levels_across_", condition_column),
      assay
    )
  }
}


# Plot Functions
# Function to add quantile lines and annotations to a violin plot
add_quantile_lines <- function(vln_plot, quantiles, colors = c("blue", "green", "red")) {
  vln_plot <- vln_plot +
    geom_hline(
      yintercept = quantiles[1],
      linetype = "dashed",
      color = colors[1]
    ) +
    geom_hline(
      yintercept = quantiles[2],
      linetype = "dashed",
      color = colors[2]
    ) +
    geom_hline(
      yintercept = quantiles[3],
      linetype = "dashed",
      color = colors[3]
    ) +
    annotate(
      "text",
      x = Inf, y = quantiles[1], label = paste0("90th: ", quantiles[1]), hjust = 1.1, color = colors[1]
    ) +
    annotate(
      "text",
      x = Inf, y = quantiles[2], label = paste0("95th: ", quantiles[2]), hjust = 1.1, color = colors[2]
    ) +
    annotate(
      "text",
      x = Inf, y = quantiles[3], label = paste0("99th: ", quantiles[3]), hjust = 1.1, color = colors[3]
    )

  return(vln_plot)
}

seurat_qc_plot <- function(seurat_obj, ncol, seurat_out_dir, prefix_file) {
  metadata_qc <- seurat_obj@meta.data
  quantile_list <- quantiles_for_qc(seurat_obj)

  nFeature_quantiles <- quantile_list[[1]]
  nCount_quantiles <- quantile_list[[2]]
  mt_quantiles <- quantile_list[[3]]
  rb_quantiles <- quantile_list[[4]]

  vp_ncount_rna <- VlnPlot(
    seurat_obj,
    features = "nCount_RNA",
    layer = "counts",
    group.by = "orig.ident",
    raster = FALSE,
    alpha = 0.2
  )
  meta_dist <- metadata_qc %>%
    ggplot2::ggplot(
      aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)
    ) +
    ggplot2::geom_density(alpha = 0.2) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_log10()

  vp <- Seurat::VlnPlot(
    seurat_obj,
    features = c(
      "nFeature_RNA", "nCount_RNA",
      "percent.mt", "percent.rb"
    ),
    group.by = "orig.ident",
    ncol = ncol,
    pt.size = 0.1
  ) & theme(plot.title = element_text(size = 16))

  # Add quantile lines to each plot using the add_quantile_lines function
  vp[[1]] <- add_quantile_lines(vp[[1]], nFeature_quantiles)
  vp[[2]] <- add_quantile_lines(vp[[2]], nCount_quantiles)
  vp[[3]] <- add_quantile_lines(vp[[3]], mt_quantiles)
  vp[[4]] <- add_quantile_lines(vp[[4]], rb_quantiles)

  fs_1 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt",
    group.by = "orig.ident"
  )
  fs_2 <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    group.by = "orig.ident"
  )
  fs_3 <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "percent.rb",
    group.by = "orig.ident"
  )
  fs_4 <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "percent.rb",
    feature2 = "percent.mt",
    group.by = "orig.ident"
  )
  fs_all <- cowplot::plot_grid(
    fs_1, fs_2, fs_3, fs_4,
    labels = c("A", "B", "C", "D"),
    ncol = 2,
    align = c("h", "v"),
    label_size = 20
  )

  pdf(
    file = paste0(seurat_out_dir, "/", prefix_file, "_qc_plots.pdf"),
    height = 14.4, width = 27.3
  )
  print(vp)
  print(fs_all)
  print(vp_ncount_rna)
  print(meta_dist)
  invisible(dev.off())
}


seurat_plot_pca <- function(
    seurat_obj,
    seurat_out_dir,
    prefix,
    layer_column = NULL,
    condition_column = NULL) {
  pdf(
    file = paste0(seurat_out_dir, "/", prefix, "_pca_plots.pdf"),
    height = 14.4,
    width = 27.3
  )
  print(
    Seurat::VizDimLoadings(
      seurat_obj,
      dims = 1:2,
      reduction = "pca"
    ) &
      theme(
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 8, face = "bold")
      )
  )
  print(
    Seurat::DimPlot(
      seurat_obj,
      reduction = "pca",
      group.by = layer_column,
      split.by = condition_column
    )
  )
  print(
    Seurat::DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)
  )
  print(
    Seurat::DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)
  )
  print(Seurat::ElbowPlot(seurat_obj))
  invisible(dev.off())
}

seurat_plot_umap <- function(
    seurat_obj,
    seurat_out_dir,
    prefix,
    group_by = NULL,
    integration_method = NULL,
    split_by = NULL) {
  if (is.null(integration_method)) {
    reduction <- "pca"
    reduction_name <- "umap"
  } else {
    reduction <- method_dict[[integration_method]][2]
    reduction_name <- paste0("umap.", reduction)
  }

  pdf(
    file = paste0(seurat_out_dir, "/", prefix, "_umap_plots.pdf"),
    height = 14.4,
    width = 27.3
  )
  print(
    DimPlot(
      seurat_obj,
      reduction = reduction_name,
      group.by = group_by,
      label = TRUE,
    )
  )

  if (!is.null(split_by)) {
    print(
      DimPlot(
        seurat_obj,
        reduction = reduction_name,
        group.by = group_by,
        split.by = split_by,
        label = TRUE,
      )
    )
  }

  invisible(dev.off())
}


plot_markers <- function(
    seurat_obj,
    assay,
    seurat_out_dir,
    prefix,
    markers,
    group_by,
    split_by) {
  seurat_obj@meta.data[[split_by]] <- factor(seurat_obj@meta.data[[split_by]])

  top1 <- unique(markers[[2]]$gene)
  top6 <- unique(markers[[3]]$gene)
  top10 <- unique(markers[[4]]$gene)
  if (is.null(split_by)) {
    split_by <- NULL
    filename <- paste0(
      seurat_out_dir, "/", prefix, "_", group_by, "_all_markers.pdf"
    )
  } else {
    filename <- paste0(
      seurat_out_dir, "/", prefix, "_", group_by, "_", split_by, "_all_markers.pdf"
    )
  }

  l_top1 <- length(top1)
  max_genes_per_page <- ifelse(is.null(split_by), 5, 2)

  pdf(file = filename, height = 14.4, width = 27.3)

  genes_plot <- 1
  while (genes_plot <= l_top1) {
    indices_to_use <- genes_plot:min(genes_plot + max_genes_per_page, l_top1)

    genes_to_plot <- top1[indices_to_use]
    print(
      FeaturePlot(
        seurat_obj,
        features = genes_to_plot,
        split.by = split_by,
        ncol = 2
      )
    )
    genes_plot <- genes_plot + max_genes_per_page + 1
  }

  print(
    DoHeatmap(
      seurat_obj,
      features = top10,
      assay = assay,
      cells = 1:500,
      group.by = group_by,
      combine = TRUE,
      raster = TRUE
    ) + NoLegend()
  )

  print(
    DotPlot(
      seurat_obj,
      features = rev(as.character(top6)),
      group.by = group_by,
      split.by = split_by
    ) +
      coord_flip() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
  )
  invisible(dev.off())
}

get_idents_column_name <- function(seurat_obj) {
  # Extract the active Idents
  idents <- as.character(Idents(seurat_obj))

  # Loop over the meta.data columns
  for (column_name in colnames(seurat_obj@meta.data)) {
    # Check if the meta.data column matches the active Idents
    if (all(seurat_obj[[column_name]] == idents)) {
      return(column_name) # Return the matching column name
    }
  }

  # If no match is found, return a message
  return("No matching column found for Idents in meta.data")
}

# Example usage:
# matching_column <- get_idents_column_name(seurat_obj)
# print(matching_column)

# plot de genes
plot_de_genes <- function(
    seurat_obj,
    de_table,
    ident_column_to_use,
    condition_column,
    seurat_out_dir,
    prefix,
    assay = "RNA") {
  original_ident_column <- get_idents_column_name(seurat_obj)

  SeuratObject::Idents(seurat_obj) <- ident_column_to_use

  within_cluster_de <- "Celltype" %in% names(de_table)
  if (within_cluster_de) {
    contrast_list <- unique(de_table[["Celltype"]])
    filter_column <- "Celltype"
  } else {
    contrast_list <- unique(de_table[["Comparison"]])
    filter_column <- "Comparison"
  }


  violin_plot_de_genes <- list()
  for (contrast in contrast_list) {
    de_genes <- de_table %>%
      dplyr::filter(!!sym(filter_column) == contrast) %>%
      dplyr::filter(p_val_adj < 0.05) %>%
      dplyr::arrange(desc(avg_log2FC)) %>%
      dplyr::pull(Gene) %>%
      head(10)

    if (!(length(de_genes) >= 1)) {
      print(de_genes)
      next
    }

    idents_to_use <- ifelse(
      within_cluster_de,
      contrast,
      as.vector(stringr::str_split(contrast, "_Vs_"))
    )[[1]]

    ncolumns <- ifelse(length(de_genes) < 5, length(de_genes), 5)
    vln_plots <- VlnPlot(
      seurat_obj,
      features = de_genes,
      idents = idents_to_use,
      group.by = ident_column_to_use,
      split.by = condition_column,
      raster = TRUE,
      ncol = ncolumns,
      combine = FALSE
    )

    vln_plots <- lapply(vln_plots, function(p) {
      p + theme(
        axis.title.x = element_blank(),
        legend.position = "none"
      )
    })

    # Combine the modified plots using patchwork
    combined_plot <- wrap_plots(vln_plots, ncol = ncolumns) +
      plot_layout(guides = "collect") + # Collect all legends into one
      theme(legend.position = "bottom") # Place the unified legend at the bottom

    violin_plot_de_genes[[as.character(contrast)]] <- combined_plot
  }

  SeuratObject::Idents(seurat_obj) <- original_ident_column

  file2save <- paste0(seurat_out_dir, "/de_", prefix, ".pdf")
  pdf(file = file2save, height = 14.4, width = 27.3)
  for (contrast in contrast_list) {
    print(violin_plot_de_genes[[as.character(contrast)]])
  }

  invisible(dev.off())
}
