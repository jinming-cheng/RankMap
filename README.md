
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RankMap <img src="man/figures/RankMap.png" align="right" alt="" width="120" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/jinming-cheng/RankMap/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/jinming-cheng/RankMap/actions)
[![Codecov test
coverage](https://codecov.io/gh/jinming-cheng/RankMap/graph/badge.svg)](https://app.codecov.io/gh/jinming-cheng/RankMap)

<!-- badges: end -->

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

Quick start guide can be found
[here](https://jinming-cheng.github.io/RankMap/index.html).

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

## Citation

Please cite this article if you use RankMap:

    To cite RankMap in publications, please use:

      Cheng J, Li S, Kim S, Ang C, Chew S, Chow P, Liu N (2026). "RankMap:
      Rank-based reference mapping for fast and robust cell type annotation
      in spatial and single-cell transcriptomics." _bioRxiv_.
      doi:10.64898/2026.03.01.708931
      <https://doi.org/10.64898/2026.03.01.708931>.

    A BibTeX entry for LaTeX users is

      @Article{,
        title = {RankMap: Rank-based reference mapping for fast and robust cell type annotation in spatial and single-cell transcriptomics},
        author = {Jinming Cheng and Shengdi Li and Serim Kim and Chow Hiang Ang and Sin Chi Chew and Pierce Kah Hoe Chow and Nan Liu},
        journal = {bioRxiv},
        year = {2026},
        doi = {10.64898/2026.03.01.708931},
      }
