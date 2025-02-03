#!/opt/sw/bioinfo-tools/sources/anaconda3/envs/d_seurat510/bin/Rscript

# Function to load libraries silently
load_library <- function(x) {
  invisible(
    suppressWarnings(
      suppressPackageStartupMessages(
        library(
          x,
          character.only = TRUE,
          quietly = TRUE,
          warn.conflicts = FALSE
        )
      )
    )
  )
}

# List of required libraries
libraries2use <- c(
  "hdf5r", "leidenAlg", "igraph", "patchwork", "ggplot2",
  "Seurat", "optparse", "dplyr", "remotes", "celldex", "SingleR",
  "writexl", "seuratter", "Seurat", "SeuratObject", "harmony",
  "seuratHelper", "this.path", "future", "comprehenr", "DESeq2",
  "yaml"
)


# Load all libraries
invisible(sapply(libraries2use, load_library))

# get directory for the present file
script_path <- this.path::this.path()

source(paste0(dirname(script_path), "/useful.R", collapse = ""))
source(paste0(dirname(script_path), "/generalized_vdj.R", collapse = ""))
# Set options for command-line arguments
option_list <- list(
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "count matrix data h5 file name  [default= %default]",
    dest = "count_matrix_file"
  ),
  make_option(
    c("-a", "--aggregate-csv"),
    type = "character",
    help = "aggregate csv file [default= %default]- file with sample_id used to aggregate using cellranger. Order must be same as in cellranger aggregate. Additional columns with information about donor, condition etc should be supplied here",
    dest = "aggregate_csv_file"
  ),
  make_option(
    c("-d", "--directory"),
    type = "character",
    help = "count matrix data directory name [default= %default]",
    dest = "count_matrix_dir"
  ),
  make_option(
    c("-o", "--out-directory"),
    type = "character",
    default = "seurat_out",
    help = "output directory [default= %default]",
    dest = "outdir"
  ),
  make_option(
    c("--min-cells"),
    type = "integer",
    default = as.integer(3),
    help = "minimum cells [default= %default]",
    dest = "min.cells"
  ),
  make_option(
    c("--min-features"),
    type = "integer",
    default = as.integer(100),
    help = "minimum features [default= %default]",
    dest = "min.features"
  ),
  make_option(
    c("--max-features"),
    type = "integer", default = as.integer(3000),
    help = "maximum features [default= %default]",
    dest = "max.features"
  ),
  make_option(
    c("--percent-mt"),
    type = "double",
    # default = as.double(100),
    help = "threshold percent mitochondrial [default= %default]
    - default filtering is done using 95th quantile",
    dest = "cutoff.percent.mt"
  ),
  make_option(
    c("--percent-rb"),
    type = "double",
    # default = as.double(100),
    help = "threshold percent ribosomal [default= %default]
    - default filtering is done using 95th quantile",
    dest = "cutoff.percent.rb"
  ),
  make_option(
    c("--project"),
    type = "character",
    default = "singleCell",
    help = "output file name [default= %default]",
    dest = "project"
  ),
  make_option(
    c("--vdj-t"),
    type = "character",
    help = "V(D)J-T annotations [default= %default]",
    dest = "vdj_t"
  ),
  make_option(
    c("--vdj-b"),
    type = "character",
    help = "V(D)J-B annotations [default= %default]",
    dest = "vdj_b"
  ),
  make_option(
    c("--layer-column"),
    type = "character",
    help = "describes experimental batches, donors, or conditions [default= %default]"
  ),
  make_option(
    c("--condition-column"),
    type = "character",
    help = "main condition for comparision  [default= %default](can be same as batch variable)"
  ),
  make_option(
    c("--integration-method"),
    type = "character",
    help = "integration method  [default= %default](CCAIntegration, RPCAIntegration, HarmonyIntegration, FastMNNIntegration, scVIIntegration)",
    dest = "integration_method"
  ),
  make_option(
    c("--enable-SCTransform"),
    type = "logical",
    default = TRUE,
    help = "sctransform normalization [default= %default]",
    dest = "enable_sct"
  ),
  make_option(
    c("--perform-DE"),
    type = "logical",
    default = FALSE,
    help = "differential expression analysis [default= %default]",
    dest = "perform_de"
  ),
  make_option(
    c("--species"),
    type = "character",
    default = "human",
    help = "annotation for species [default= %default]",
    dest = "species"
  ),
  make_option(
    c("--num-cores"),
    type = "integer",
    default = as.integer(4L),
    help = "number of cores to use [default= %default]",
    dest = "num_cores"
  ),
  make_option(
    c("--num-threads"),
    type = "integer",
    default = as.integer(16L),
    help = "number of cores to use [default= %default]",
    dest = "num_threads"
  ),
  make_option(
    c("--memory-usage"),
    type = "double",
    default = as.double(16.0),
    help = "memory available to use in [default= %default Gb]",
    dest = "memory_usage"
  )
)

# Parse command-line arguments
opt_parser <- OptionParser(
  option_list = option_list,
  description = "This script processes single-cell RNA sequencing data. Performs quality control, filtering, normalization, batch-correction(optional), clustering and annotation. Optionally it also does differential expression analysis. It integrates various Seurat functions and provides pertinent figures, and tables for an exhaustive investigation."
)
opt <- parse_args(opt_parser)

# Define user inputs
project_name <- opt$project
seurat_out_dir <- opt$outdir
min_cells <- opt$min.cells
min_features <- opt$min.features
max_features <- opt$max.features
percent_mt <- opt$cutoff.percent.mt
percent_rb <- opt$cutoff.percent.rb
data_dir <- opt$count_matrix_dir
data_file <- opt$count_matrix_file
aggr_csv_file <- opt$aggregate_csv_file
tcr_file <- opt$vdj_t
bcr_file <- opt$vdj_b
layer_column <- opt$layer_column
condition_column <- opt$condition_column
integration_method <- opt$integration_method
enable_sct <- opt$enable_sct
perform_de <- opt$perform_de
species <- opt$species
cores <- opt$num_cores
threads <- opt$num_threads
memory <- opt$memory_usage


# resources
RhpcBLASctl::blas_set_num_threads(threads)
RhpcBLASctl::omp_set_num_threads(threads)

memory_to_use <- memory * 1024^3

options(future.globals.maxSize = memory_to_use)
options(mc.cores = cores)

invisible(plan())

# Main Seurat analysis function
seurat_analysis <- function(
    data_dir = NULL,
    data_file = NULL,
    project_name = "singleCell",
    seurat_out_dir = "seurat_out",
    min_cells = 3L,
    min_features = 100L,
    max_features = 3000L,
    percent_mt = NULL,
    percent_rb = NULL,
    aggr_csv_file = NULL,
    tcr_file = NULL,
    bcr_file = NULL,
    layer_column = NULL,
    condition_column = NULL,
    integration_method = NULL,
    enable_sct = TRUE,
    perform_de = FALSE,
    species = "human") {
  start_time <- Sys.time()
  if (!dir.exists(seurat_out_dir)) {
    dir.create(seurat_out_dir)
  }

  # Load dataset
  data_10_x <- load_data_10X(data_dir, data_file)
  print("loaded data_10X successfully")
  # Load metadata if applicable
  sample_identity <- metadata(data_10_x, aggr_csv_file)

  # Create Seurat object with multiple layers
  seurat_obj <- create_seurat_object(
    data_10_x,
    sample_identity,
    project_name
  )

  seurat_obj <- seurat_set_orig_ident(seurat_obj, sample_identity)

  # Add V(D)J-T and B annotations if applicable
  seurat_obj <- add_clonotype(tcr_file, seurat_obj, "t")
  seurat_obj <- add_clonotype(bcr_file, seurat_obj, "b")


  # Quality control
  seurat_obj <- seurat_qc(seurat_obj)

  seurat_qc_plot(seurat_obj, ncol = 4, seurat_out_dir, "pre_filter")

  # Save Raw Seurat object
  file2save <- paste0(seurat_out_dir, "/", project_name, "_raw_seurat.rds")
  saveRDS(
    seurat_obj,
    file = file2save
  )

  seurat_obj <- seurat_filtering(
    seurat_obj,
    min_cells,
    min_features,
    max_features,
    percent_mt,
    percent_rb
  )
  seurat_qc_plot(seurat_obj, ncol = 4, seurat_out_dir, "post_filter")
  # split if  datasets ran across experimental batches, donors, or condition
  seurat_obj <- split_seurat_by_layer(
    seurat_obj,
    assay = "RNA",
    layer_column = layer_column
  )

  message("Number of parallel workers: ", nbrOfWorkers())

  seurat_obj <- seurat_normalization(seurat_obj, sct = enable_sct)
  seurat_obj <- seurat_dim_reduction_linear(seurat_obj)
  seurat_plot_pca(
    seurat_obj,
    seurat_out_dir,
    project_name,
    layer_column = layer_column,
    condition_column = condition_column
  )

  if (!is.null(layer_column)) {
    seurat_obj <- seurat_dim_reduction_nonlinear(
      seurat_obj,
      integration_method = NULL
    )
    seurat_plot_umap(
      seurat_obj,
      seurat_out_dir,
      paste0(project_name, "_no_integration"),
      group_by = layer_column,
      integration_method = NULL,
      split_by = condition_column
    )
  }

  # Integrative analysis using the specified method
  seurat_obj <- integrate_seurat_layers(
    seurat_obj = seurat_obj,
    integration_method = integration_method,
    enable_sct = enable_sct
  )

  if (!is.null(integration_method)) {
    seurat_obj <- seurat_dim_reduction_nonlinear(
      seurat_obj,
      integration_method = integration_method
    )

    seurat_plot_umap(
      seurat_obj,
      seurat_out_dir,
      paste0(project_name, "_", integration_method, "_integrated"),
      group_by = layer_column,
      integration_method = integration_method,
      split_by = condition_column
    )
  }

  assay <- SeuratObject::DefaultAssay(seurat_obj)

  seurat_obj <- join_seurat_layers(
    seurat_obj,
    assay = "RNA",
    layer_column = layer_column
  )

  SeuratObject::DefaultAssay(seurat_obj) <- assay
  seurat_obj <- seurat_clustering(seurat_obj, integration_method)
  seurat_obj <- seurat_dim_reduction_nonlinear(seurat_obj)


  seurat_plot_umap(
    seurat_obj,
    seurat_out_dir,
    prefix = paste0(
      project_name,
      ifelse(is.null(layer_column), "", paste0("_", layer_column, "_integrated")), "_seurat_clusters_final"
    ),
    group_by = NULL,
    integration_method = integration_method,
    split_by = NULL
  )

  seurat_obj <- seurat_annotation(seurat_obj, species = species)

  if ("seurat_clusters" %in% names(seurat_obj@meta.data)) {
    group_by <- c("seurat_clusters", "SingleR.labels")
  } else {
    group_by <- c("SingleR.labels")
  }

  seurat_plot_umap(
    seurat_obj,
    seurat_out_dir,
    prefix = paste0(project_name, "_annotated"),
    group_by = group_by,
    integration_method = integration_method,
    split_by = condition_column
  )

  if (perform_de) {
    # revert to "RNA" assay for DE analysis
    SeuratObject::DefaultAssay(seurat_obj) <- "RNA"

    if (sum(c("data", "scale.data") %in% SeuratObject::Layers(seurat_obj@assays$RNA)) == 2) {
      cat("The \"RNA\" assay is already normalized.\n")
    } else {
      cat("The RNA assay is not normalized.\n")
      cat("Normalizing RNA assay ...\n")
      seurat_obj <- seurat_normalization(seurat_obj, sct = FALSE)
    }

    for (group_column in group_by) {
      cat(paste0("processing ", group_column, " for DE\n"))

      seurat_markers <- seurat_all_markers(
        seurat_obj,
        assay = "RNA",
        cluster_column = group_column,
        prefix = project_name,
        seurat_out_dir = seurat_out_dir
      )

      plot_markers(
        seurat_obj,
        "RNA",
        seurat_out_dir,
        project_name,
        seurat_markers,
        group_column,
        condition_column
      )

      seurat_de(
        seurat_obj,
        annotation_column = group_column,
        assay = "RNA",
        seurat_out_dir
      )

      seurat_de_across_condition_same_cell_type(
        seurat_obj,
        annotation_column = group_column,
        condition_column = condition_column,
        assay = "RNA",
        seurat_out_dir
      )

      groups_to_use <- if ("donor" %in% names(seurat_obj@meta.data)) {
        c(condition_column, "donor", group_column)
      } else {
        c(condition_column, group_column)
      }

      seurat_de_pseudo_bulk_across_condition(
        seurat_obj,
        annotation_column = group_column,
        condition_column = condition_column,
        group_by = groups_to_use,
        assay = "RNA",
        seurat_out_dir
      )
    }
  }

  # Save Seurat object
  file2save <- paste0(
    seurat_out_dir, "/", project_name, "_analysed_seurat.rds"
  )

  saveRDS(
    seurat_obj,
    file = file2save
  )

  end_time <- Sys.time()
  print("Seurat analysis complete!!")
  print(end_time - start_time)

  sink(paste0(seurat_out_dir, "/session_info.txt"))
  sessionInfo()
  sink()
}


# data_dir <- "count/filtered_feature_bc_matrix"
# # data_dir <- "count/sample_filtered_feature_bc_matrix"
# data_file <- NULL
# project_name <- "test"
# seurat_out_dir <- "seurat_out"
# min_cells <- 3
# min_features <- 100
# max_features <- 3000
# percent_mt <- NULL
# percent_rb <- NULL
# aggr_csv_file <- "aggregation.csv"
# tcr_file <- "vdj_t/filtered_contig_annotations.csv"
# bcr_file <- "vdj_b/filtered_contig_annotations.csv"
# layer_column <- "donor"
# condition_column <- "health_status"
# integration_method <- "RPCAIntegration"
# # CCAIntegration, RPCAIntegration, HarmonyIntegration,
# # FastMNNIntegration, scVIIntegration
# enable_sct <- TRUE
# perform_de <- FALSE
# species <- "human"

# Run Seurat analysis
if ((!is.null(data_dir)) | (!is.null(data_file))) {
  seurat_analysis(
    data_dir,
    data_file,
    project_name,
    seurat_out_dir,
    min_cells,
    min_features,
    max_features,
    percent_mt,
    percent_rb,
    aggr_csv_file,
    tcr_file,
    bcr_file,
    layer_column,
    condition_column,
    integration_method,
    enable_sct,
    perform_de,
    species
  )
}
