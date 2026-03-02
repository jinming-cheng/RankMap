# RankMap ![](reference/figures/RankMap.png)

RankMap is a fast, robust, and scalable method for reference-based cell
type annotation of single-cell and spatial transcriptomics data. It
transforms gene expression matrices into sparse ranked representations
and applies multinomial regression via the glmnet framework to predict
cell types. This rank-based strategy enhances robustness to batch
effects, platform differences, and partial gene coverage, making RankMap
especially suitable for technologies like Xenium and MERFISH. The
package supports flexible preprocessing, fast model training and
prediction, and is compatible with Seurat, SingleCellExperiment, and
SpatialExperiment objects. RankMap achieves competitive accuracy with
significantly lower runtime than existing methods such as SingleR,
Azimuth, and RCTD.

## Installation

The *RankMap* package can be installed from GitHub by using:

``` r
devtools::install_github("jinming-cheng/RankMap")
```

## Quick Start

An example of predicting cell types for Xenium spatial data using a
well-annotated single-cell dataset as reference

``` r
# load data
seu_sc <- readRDS(system.file("extdata", "seu_sc.rds",
    package = "RankMap"
))

seu_xen <- readRDS(system.file("extdata", "seu_xen.rds",
    package = "RankMap"
))

# predict cell types, if unsure about k (default 20), set k to 100.
pred_df <- RankMap(
    ref_data = seu_sc,
    ref_labels = seu_sc$cell_type,
    new_data = seu_xen
    k = 100
)
head(pred_df)
```
