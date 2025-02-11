# scQCAD

[![Repository](https://img.shields.io/badge/repository-github-blue?color=%6495ED)](https://github.com/BioDebojyoti/scQCAD)
[![Author](https://img.shields.io/badge/author-Debojyoti%20Das-purple?color=%9FE2BF)](https://github.com/BioDebojyoti/BioDebojyoti)

This Seurat wrapper takes aggregated cellranger output, and performs:

- Quality control,
- Batch correction (optional),
- Clustering
- Cell-type annotation using [SingleR](https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html)
- Differential expression (DE) analysis (optional)  
  using the Seurat package in R. It takes as input a sample sheet specifying metadata and directory paths for each sample and generates a range of outputs for QC, visualization, and downstream analysis.

## Input Requirements

The wrapper will need the [Cell Ranger v8.0.1](https://github.com/10XGenomics/cellranger/releases/tag/cellranger-8.0.1) output folder path or the h5 file following **aggr** step. It will also need a sample sheet (aggregate_csv_file) containing the meta-data for the different samples.
The aggregate_csv_file should be a CSV file with the following columns (it should follow the same order as used in cellranger aggr step).

| sample_id | sample_outs        | donor | origin     | group   |
| :-------- | :----------------- | :---- | :--------- | :------ |
| S2        | "/path/to/outs/S2" | S2    | S2_Disease | Disease |
| S9        | "/path/to/outs/S9" | S9    | S9_Normal  | Normal  |
| S1        | "/path/to/outs/S1" | S1    | S1_Disease | Disease |
| S4        | "/path/to/outs/S4" | S4    | S4_Disease | Disease |

## Output Files

scQCAD will produce figures in pdf format, and tables as CSV files. It will also produce two R objects in output. The R objects can be used later for additional or (re-analysis) analysisis. Here is a brief description of the output files content:

| filename                                                                                                                                                                                                                                                                | description                                                                                                                                                                        |
| :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Quality control (QC)**                                                                                                                                                                                                                                                |                                                                                                                                                                                    |
| <span style="color:red">pre_filter</span>\_qc_plot.pdf                                                                                                                                                                                                                  | pre-filtering QC metrics                                                                                                                                                           |
| <span style="color:green">post_filter</span>\_qc_plot.pdf                                                                                                                                                                                                               | post-filtering QC metrics                                                                                                                                                          |
| **Principal Component Analysis (PCA)**                                                                                                                                                                                                                                  |                                                                                                                                                                                    |
| **\<project\>**\_<span style="color:green">pca</span>\_plots.pdf                                                                                                                                                                                                        | PCA plots                                                                                                                                                                          |
| **Uniform Manifold Approximation and Projection (UMAP)**                                                                                                                                                                                                                |                                                                                                                                                                                    |
| **\<project\>**\_<span style="color:orange">no_integration</span>\_umap_plots.pdf                                                                                                                                                                                       | UMAP (conditonal on presence of batch variable) before integration                                                                                                                 |
| **Integration (UMAP)**                                                                                                                                                                                                                                                  |                                                                                                                                                                                    |
| **\<project\>**\_**\<integration_method\>**\_<span style="color:green">integrated</span>\_umap_plots.pdf                                                                                                                                                                | UMAP (conditonal on presence of batch variable) after integration                                                                                                                  |
| **Clustering (UMAP)**                                                                                                                                                                                                                                                   |                                                                                                                                                                                    |
| **\<project\>**\_seurat_clusters\_<span style="color:green">final</span>\_umap_plots.pdf <br>OR</br> **\<project\>**\_**\<layer_column\>**\_<span style="color:green">integrated</span>\_seurat_clusters\_<span style="color:green">final</span>\_umap_plots.pdf</span> | umap plot(s) of seurat clusters (conditonal on presence of batch variable **\<layer_column\>** after integration)                                                                  |
| **Annotation (UMAP)**                                                                                                                                                                                                                                                   |                                                                                                                                                                                    |
| **\<project\>**\_<span style="color:green">annotated</span>\_umap_plots.pdf                                                                                                                                                                                             | umap plot(s) of **singleR\*** annotated clusters (split by **\<condition_column\>** if present)                                                                                    |
| **Marker Identification**                                                                                                                                                                                                                                               |                                                                                                                                                                                    |
| **\<project\>**\_**\<cluster_name\>**\_<span style="color:green">all_markers</span>.pdf <br>OR</br> **\<project\>**\_**\<cluster_name\>**\_**\<condition_column\>**\_<span style="color:green">all_markers</span>.pdf                                                   | FeaturePlot of top1 markers, Heatmap of top10 markers, and DotPlot of top6 markers of **\<cluster_name\>**. FeaturePlot and DotPlots split by **\<condition_column\>** if present. |

**\<project\>** default name is "singleCell"\
**\<integration_method\>** no default integration-method\
**\<layer_column\>** no default batch-variable\
**\<condition_column\>** no default group-variable\
**singleR\*** currently only supports species "human", and "mouse".\
**\<cluster_name\>** "seurat_clusters" and "singleR.labels".

Currently, the following references are used for cell type annotation

```r
# species = "human"
celldex::HumanPrimaryCellAtlasData()
# species = "mouse"
celldex::MouseRNAseqData()
```

The wrapper run saves two R objects

```r
"/path/to/<seurat_out_dir>/<project>_raw_seurat.rds"
"/path/to/<seurat_out_dir>/<project>_analysed_seurat.rds"
```

The former can be used to re-analyse without having to create the seurat object afresh. The latter object is meant to help perform additional analysis should there be a need for it. Please follow the steps described below to load the object of interest.

```r
library(Seurat)
# for re-analysis
seurat_obj <- readRDS("path/to/<seurat_out_dir>/<project>_raw_seurat.rds")
# for additional analysis
seurat_obj <- readRDS("path/to/<seurat_out_dir>/<project>_analysed_seurat.rds")
```

While the wrapper is expected to be used by command-line invocation it can also be used as a function by first making it available as shown below:

## Make the functions available

```r
source("scQCAD.R")
```

## Run the analysis

```r
# Run Seurat analysis
seurat_analysis(
data_dir = "count/filtered_feature_bc_matrix",
data_file = NULL,
project_name = "test",
seurat_out_dir = "seurat_out",
min_cells = 3,
min_features = 100,
max_features = 3000,
percent_mt = NULL,
percent_rb = NULL,
aggr_csv_file = "aggregation.csv",
tcr_file = "vdj_t/filtered_contig_annotations.csv",
bcr_file = "vdj_b/filtered_contig_annotations.csv",
layer_column = "donor",
condition_column = "health_status",
integration_method = "RPCAIntegration",
enable_sct = TRUE,
perform_de = FALSE,
species = "human"
)
```

## CLI options

```bash
$ Rscript scQCAD.R --help
Usage: scQCAD.R [options]
This script processes single-cell RNA sequencing data. Performs quality control,
filtering, normalization, batch-correction(optional), clustering and annotation.
Optionally it also does differential expression analysis. It integrates various
Seurat functions and provides pertinent figures, and tables for an exhaustive
investigation.

Options:
	-f FILE, --file=FILE
		count matrix data h5 file name  [default= NULL]

	-a AGGREGATE-CSV, --aggregate-csv=AGGREGATE-CSV
		aggregate csv file [default= NULL] - file with sample_id used to aggregate
      using cellranger. Order must be same as in cellranger aggregate. Additional
      columns with information about donor, condition etc should be supplied here

	-d DIRECTORY, --directory=DIRECTORY
		count matrix data directory name [default= NULL]

	-o OUT-DIRECTORY, --out-directory=OUT-DIRECTORY
		output directory [default= seurat_out]

	--min-cells=MIN-CELLS
		minimum cells [default= 3]

	--min-features=MIN-FEATURES
		minimum features [default= 100]

	--max-features=MAX-FEATURES
		maximum features [default= 3000]

	--percent-mt=PERCENT-MT
		threshold percent mitochondrial [default= NULL] - default filtering is done
      using 95th quantile

	--percent-rb=PERCENT-RB
		threshold percent ribosomal [default= NULL] - default filtering is done using
      95th quantile

	--project=PROJECT
		output file name [default= singleCell]

	--vdj-t=VDJ-T
		V(D)J-T annotations [default= NULL]

	--vdj-b=VDJ-B
		V(D)J-B annotations [default= NULL]

	--layer-column=LAYER-COLUMN
		describes experimental batches, donors, or conditions [default= NULL]

	--condition-column=CONDITION-COLUMN
		main condition for comparision  [default= NULL]

	--integration-method=INTEGRATION-METHOD
		integration method  [default= NULL] - (CCAIntegration, RPCAIntegration,
      HarmonyIntegration, FastMNNIntegration, scVIIntegration)

	--enable-SCTransform=ENABLE-SCTRANSFORM
		sctransform normalization [default= TRUE]

	--perform-DE=PERFORM-DE
		differential expression analysis [default= FALSE]

	--species=SPECIES
		annotation for species [default= human]

	-h, --help
		Show this help message and exit
```

## Differential Expression (DE) Analysis (optional)

DE analysis is optional and is performed if \-\-perform\-DE is enabled. Three different scenarios are explored:

- DE between each pair of clusters (seurat_clusters and SingleR.label)
- DE analysis within the same cell type (cluster type) across conditions
- DE analysis with pseudo bulking
