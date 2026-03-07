# Getting Started with RankMap

## Introduction

**RankMap** is an R package for fast, robust, and scalable
reference-based cell type annotation in single-cell and spatial
transcriptomics data. It works by transforming gene expression matrices
into sparse ranked representations and training a multinomial logistic
regression model using the `glmnet` framework. This rank-based approach
improves robustness to batch effects, platform differences, and partial
gene coverage—especially beneficial for technologies such as Xenium and
MERFISH.

**RankMap** supports commonly used data structures including `Seurat`,
`SingleCellExperiment`, and `SpatialExperiment`. The workflow includes
flexible preprocessing steps such as top-K gene masking, binning,
expression weighting, and scaling, followed by efficient model training
and rapid prediction.

Compared to existing tools such as **SingleR**, **RCTD** (via
**spacexr**), and **Azimuth**, **RankMap** achieves comparable or
superior accuracy with significantly faster runtime, making it
particularly well suited for high-throughput applications on large
datasets.

This vignette provides a quick-start guide to using **RankMap** for cell
type prediction.

## Installation

Install **RankMap** from GitHub:

``` r
devtools::install_github("jinming-cheng/RankMap")
```

## Quick Start (Seurat Objects)

### Load Data

``` r
library(RankMap)
```

Load example single-cell RNA-seq dataset (17,597 genes x 150 cells):

``` r
seu_sc <- readRDS(system.file("extdata", "seu_sc.rds", package = "RankMap"))
seu_sc
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 'SeuratObject' was built under R 4.4.0 but the current version is
#> 4.4.2; it is recomended that you reinstall 'SeuratObject' as the ABI
#> for R may have changed
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
#> An object of class Seurat 
#> 17597 features across 150 samples within 1 assay 
#> Active assay: RNA (17597 features, 0 variable features)
#>  2 layers present: counts, data
```

Load example Xenium spatial transcriptomics dataset (313 genes x 150
cells):

``` r
seu_xen <- readRDS(system.file("extdata", "seu_xen.rds", package = "RankMap"))
seu_xen
#> An object of class Seurat 
#> 313 features across 150 samples within 1 assay 
#> Active assay: RNA (313 features, 0 variable features)
#>  2 layers present: counts, data
```

### Predict Cell Types

Run cell type prediction using the
[`RankMap()`](https://github.com/jinming-cheng/RankMap/reference/RankMap.md)
function. By default, RankMap uses normalized expression from the “data”
slot. For spatial datasets with limited gene panels, a smaller `k`
(e.g., `k = 20`) is typically sufficient. For single-cell RNA-seq with
deeper coverage, larger values of `k` (e.g., 100 or 200) are generally
recommended.

``` r
pred_df <- RankMap(
    ref_data = seu_sc,
    ref_labels = seu_sc$cell_type,
    new_data = seu_xen,
    k = 20
)
```

The result is a `data.frame` containing: `cell_id`,
`predicted_cell_type` and `confidence`

``` r
head(pred_df)
#>   cell_id predicted_cell_type confidence
#> 1    3869               Tumor     0.8829
#> 2    5257               Tumor     0.9612
#> 3    6456               Basal     0.9243
#> 4    8555                  LP     0.8847
#> 5    9243               Basal     0.9911
#> 6   10303               Basal     0.9971
```

### Evaluate Performance

If ground truth labels are available, you can evaluate prediction
accuracy using:

``` r
perf <- evaluatePredictionPerformance(
    prediction_df = pred_df,
    truth = seu_xen$cell_type_SingleR
)
perf
#> $overall_accuracy
#> [1] 0.9466667
#> 
#> $per_class_accuracy
#> Basal    LP Tumor 
#>  0.96  0.90  0.98 
#> 
#> $confusion_matrix
#>        Predicted
#> True    Basal LP Tumor
#>   Basal    48  2     0
#>   LP        4 45     1
#>   Tumor     1  0    49
```

## Quick Start (SummarizedExperiment Objects)

### Prepare Data

Convert `Seurat` objects into `SingleCellExperiment` objects:

``` r
library(SingleCellExperiment)
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following object is masked from 'package:SeuratObject':
#> 
#>     intersect
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
#>     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
#>     table, tapply, union, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> 
#> Attaching package: 'IRanges'
#> The following object is masked from 'package:sp':
#> 
#>     %over%
#> Loading required package: GenomeInfoDb
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
#> 
#> Attaching package: 'SummarizedExperiment'
#> The following object is masked from 'package:SeuratObject':
#> 
#>     Assays

sce_sc <- SingleCellExperiment(
    assays = list(
        counts = GetAssayData(seu_sc, layer = "counts"),
        logcounts = GetAssayData(seu_sc, layer = "data")
    ),
    colData = seu_sc[[]] # seu_sc@meta.data
)

sce_sp <- SingleCellExperiment(
    assays = list(
        counts = GetAssayData(seu_xen, layer = "counts"),
        logcounts = GetAssayData(seu_xen, layer = "data")
    ),
    colData = seu_xen[[]] # seu_xen@meta.data
)
```

### Predict Cell Types

Run cell type prediction using the
[`RankMap()`](https://github.com/jinming-cheng/RankMap/reference/RankMap.md)
function. Set `k = 100` as a reasonable default when the optimal number
of top-ranked genes is unknown. When using `SummarizedExperiment` input,
the `logcounts` assay is used automatically.

``` r
pred_df <- RankMap(
    ref_data = sce_sc,
    ref_labels = sce_sc$cell_type,
    new_data = sce_sp,
    k = 100
)
```

### Evaluate Performance

Compare predictions with ground truth labels:

``` r
perf <- evaluatePredictionPerformance(
    prediction_df = pred_df,
    truth = sce_sp$cell_type_SingleR
)
perf
#> $overall_accuracy
#> [1] 0.98
#> 
#> $per_class_accuracy
#> Basal    LP Tumor 
#>  0.98  1.00  0.96 
#> 
#> $confusion_matrix
#>        Predicted
#> True    Basal LP Tumor
#>   Basal    49  1     0
#>   LP        0 50     0
#>   Tumor     2  0    48
```

## Session Info

``` r
sessionInfo()
#> R version 4.4.2 (2024-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0
#>  [3] Biobase_2.66.0              GenomicRanges_1.58.0       
#>  [5] GenomeInfoDb_1.42.3         IRanges_2.40.1             
#>  [7] S4Vectors_0.44.0            BiocGenerics_0.52.0        
#>  [9] MatrixGenerics_1.18.1       matrixStats_1.5.0          
#> [11] SeuratObject_5.3.0          sp_2.2-1                   
#> [13] RankMap_0.99.0              BiocStyle_2.34.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3      jsonlite_2.0.0          shape_1.4.6.1          
#>   [4] magrittr_2.0.4          spatstat.utils_3.2-1    farver_2.1.2           
#>   [7] rmarkdown_2.30          zlibbioc_1.52.0         fs_1.6.7               
#>  [10] ragg_1.5.1              vctrs_0.7.1             ROCR_1.0-12            
#>  [13] spatstat.explore_3.7-0  S4Arrays_1.6.0          htmltools_0.5.9        
#>  [16] SparseArray_1.6.2       sass_0.4.10             sctransform_0.4.3      
#>  [19] parallelly_1.46.1       KernSmooth_2.23-26      bslib_0.10.0           
#>  [22] htmlwidgets_1.6.4       desc_1.4.3              ica_1.0-3              
#>  [25] plyr_1.8.9              plotly_4.12.0           zoo_1.8-15             
#>  [28] cachem_1.1.0            igraph_2.2.2            mime_0.13              
#>  [31] lifecycle_1.0.5         iterators_1.0.14        pkgconfig_2.0.3        
#>  [34] Matrix_1.7-4            R6_2.6.1                fastmap_1.2.0          
#>  [37] GenomeInfoDbData_1.2.13 fitdistrplus_1.2-6      future_1.69.0          
#>  [40] shiny_1.13.0            digest_0.6.39           patchwork_1.3.2        
#>  [43] Seurat_5.4.0            tensor_1.5.1            RSpectra_0.16-2        
#>  [46] irlba_2.3.7             textshaping_1.0.5       progressr_0.18.0       
#>  [49] spatstat.sparse_3.1-0   httr_1.4.8              polyclip_1.10-7        
#>  [52] abind_1.4-8             compiler_4.4.2          S7_0.2.1               
#>  [55] fastDummies_1.7.5       MASS_7.3-65             DelayedArray_0.32.0    
#>  [58] tools_4.4.2             lmtest_0.9-40           otel_0.2.0             
#>  [61] httpuv_1.6.16           future.apply_1.20.2     goftest_1.2-3          
#>  [64] glue_1.8.0              nlme_3.1-168            promises_1.5.0         
#>  [67] grid_4.4.2              Rtsne_0.17              cluster_2.1.8.2        
#>  [70] reshape2_1.4.5          generics_0.1.4          gtable_0.3.6           
#>  [73] spatstat.data_3.1-9     tidyr_1.3.2             data.table_1.18.2.1    
#>  [76] XVector_0.46.0          spatstat.geom_3.7-0     RcppAnnoy_0.0.23       
#>  [79] ggrepel_0.9.6           RANN_2.6.2              foreach_1.5.2          
#>  [82] pillar_1.11.1           stringr_1.6.0           spam_2.11-3            
#>  [85] RcppHNSW_0.6.0          later_1.4.8             splines_4.4.2          
#>  [88] dplyr_1.2.0             lattice_0.22-9          survival_3.8-6         
#>  [91] deldir_2.0-4            tidyselect_1.2.1        miniUI_0.1.2           
#>  [94] pbapply_1.7-4           knitr_1.51              gridExtra_2.3          
#>  [97] bookdown_0.46           scattermore_1.2         xfun_0.56              
#> [100] UCSC.utils_1.2.0        stringi_1.8.7           lazyeval_0.2.2         
#> [103] yaml_2.3.12             evaluate_1.0.5          codetools_0.2-20       
#> [106] tibble_3.3.1            BiocManager_1.30.27     cli_3.6.5              
#> [109] uwot_0.2.4              xtable_1.8-8            reticulate_1.45.0      
#> [112] systemfonts_1.3.2       jquerylib_0.1.4         Rcpp_1.1.1             
#> [115] globals_0.19.0          spatstat.random_3.4-4   png_0.1-8              
#> [118] spatstat.univar_3.1-6   parallel_4.4.2          pkgdown_2.2.0          
#> [121] ggplot2_4.0.2           dotCall64_1.2           listenv_0.10.0         
#> [124] glmnet_4.1-10           viridisLite_0.4.3       scales_1.4.0           
#> [127] ggridges_0.5.7          crayon_1.5.3            purrr_1.2.1            
#> [130] rlang_1.1.7             cowplot_1.2.0
```
