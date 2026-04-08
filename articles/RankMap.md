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

Install RankMap from Bioconductor

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("RankMap")
```

## Quick Start (Seurat Objects)

### Load Data

``` r

library(RankMap)
library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 'SeuratObject' was built under R 4.6.0 but the current version is
#> 4.7.0; it is recomended that you reinstall 'SeuratObject' as the ABI
#> for R may have changed
#> 'SeuratObject' was built with package 'Matrix' 1.7.4 but the current
#> version is 1.7.5; it is recomended that you reinstall 'SeuratObject' as
#> the ABI for 'Matrix' may have changed
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
```

Load example single-cell RNA-seq dataset (17,597 genes x 150 cells):

``` r

seu_sc <- readRDS(system.file("extdata", "seu_sc.rds", package = "RankMap"))
seu_sc
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
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
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
#> Loading required package: Seqinfo
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
#> The following object is masked from 'package:Seurat':
#> 
#>     Assays
#> The following object is masked from 'package:SeuratObject':
#> 
#>     Assays
```

``` r

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
#> R Under development (unstable) (2026-04-05 r89793)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
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
#>  [1] SingleCellExperiment_1.33.2 SummarizedExperiment_1.41.1
#>  [3] Biobase_2.71.0              GenomicRanges_1.63.1       
#>  [5] Seqinfo_1.1.0               IRanges_2.45.0             
#>  [7] S4Vectors_0.49.1            BiocGenerics_0.57.0        
#>  [9] generics_0.1.4              MatrixGenerics_1.23.0      
#> [11] matrixStats_1.5.0           Seurat_5.4.0               
#> [13] SeuratObject_5.3.0          sp_2.2-1                   
#> [15] RankMap_0.99.1              BiocStyle_2.39.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         shape_1.4.6.1         
#>   [4] magrittr_2.0.5         spatstat.utils_3.2-2   farver_2.1.2          
#>   [7] rmarkdown_2.31         fs_2.0.1               ragg_1.5.2            
#>  [10] vctrs_0.7.2            ROCR_1.0-12            spatstat.explore_3.8-0
#>  [13] S4Arrays_1.11.1        htmltools_0.5.9        SparseArray_1.11.13   
#>  [16] sass_0.4.10            sctransform_0.4.3      parallelly_1.46.1     
#>  [19] KernSmooth_2.23-26     bslib_0.10.0           htmlwidgets_1.6.4     
#>  [22] desc_1.4.3             ica_1.0-3              plyr_1.8.9            
#>  [25] plotly_4.12.0          zoo_1.8-15             cachem_1.1.0          
#>  [28] igraph_2.2.3           mime_0.13              lifecycle_1.0.5       
#>  [31] iterators_1.0.14       pkgconfig_2.0.3        Matrix_1.7-5          
#>  [34] R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-6    
#>  [37] future_1.70.0          shiny_1.13.0           digest_0.6.39         
#>  [40] patchwork_1.3.2        tensor_1.5.1           RSpectra_0.16-2       
#>  [43] irlba_2.3.7            textshaping_1.0.5      progressr_0.19.0      
#>  [46] spatstat.sparse_3.1-0  httr_1.4.8             polyclip_1.10-7       
#>  [49] abind_1.4-8            compiler_4.7.0         S7_0.2.1              
#>  [52] fastDummies_1.7.5      MASS_7.3-65            DelayedArray_0.37.1   
#>  [55] tools_4.7.0            lmtest_0.9-40          otel_0.2.0            
#>  [58] httpuv_1.6.17          future.apply_1.20.2    goftest_1.2-3         
#>  [61] glue_1.8.0             nlme_3.1-169           promises_1.5.0        
#>  [64] grid_4.7.0             Rtsne_0.17             cluster_2.1.8.2       
#>  [67] reshape2_1.4.5         gtable_0.3.6           spatstat.data_3.1-9   
#>  [70] tidyr_1.3.2            data.table_1.18.2.1    XVector_0.51.0        
#>  [73] spatstat.geom_3.7-3    RcppAnnoy_0.0.23       ggrepel_0.9.8         
#>  [76] RANN_2.6.2             foreach_1.5.2          pillar_1.11.1         
#>  [79] stringr_1.6.0          spam_2.11-3            RcppHNSW_0.6.0        
#>  [82] later_1.4.8            splines_4.7.0          dplyr_1.2.1           
#>  [85] lattice_0.22-9         survival_3.8-6         deldir_2.0-4          
#>  [88] tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-4         
#>  [91] knitr_1.51             gridExtra_2.3          bookdown_0.46         
#>  [94] scattermore_1.2        xfun_0.57              stringi_1.8.7         
#>  [97] lazyeval_0.2.3         yaml_2.3.12            evaluate_1.0.5        
#> [100] codetools_0.2-20       tibble_3.3.1           BiocManager_1.30.27   
#> [103] cli_3.6.5              uwot_0.2.4             xtable_1.8-8          
#> [106] reticulate_1.45.0      systemfonts_1.3.2      jquerylib_0.1.4       
#> [109] Rcpp_1.1.1             globals_0.19.1         spatstat.random_3.4-5 
#> [112] png_0.1-9              spatstat.univar_3.1-7  parallel_4.7.0        
#> [115] pkgdown_2.2.0          ggplot2_4.0.2          dotCall64_1.2         
#> [118] listenv_0.10.1         glmnet_4.1-10          viridisLite_0.4.3     
#> [121] scales_1.4.0           ggridges_0.5.7         purrr_1.2.1           
#> [124] rlang_1.2.0            cowplot_1.2.0
```
