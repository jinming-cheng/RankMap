# Changelog

## RankMap 0.99.1

### Changes

- Updated `DESCRIPTION` to require `R (>= 4.6.0)`.
- Refactored repeated matrix transposition logic in
  [`trainRankModel()`](https://github.com/jinming-cheng/RankMap/reference/TrainRankModel.md)
  and
  [`predictRankModel()`](https://github.com/jinming-cheng/RankMap/reference/PredictRankModel.md)
  into an internal helper function (`.transposeMatrix`).
- Reduced use of `::` by importing selected functions (e.g., `glmnet()`,
  `cv.glmnet()`) to improve dependency checking.
- Added package-level documentation (`RankMap-package.R`).
- Updated vignette:
  - Removed GitHub installation instructions.
  - Added Bioconductor installation instructions using
    `BiocManager::install("RankMap")`.
  - Added unique labels to all code chunks.

## RankMap 0.99.0

*Initial submission to Bioconductor*

### New Features

- Fast, robust, and scalable reference-based cell type annotation using
  multinomial regression on ranked expression matrices.
- Supports both **single-cell** and **spatial transcriptomics** data.
- Compatible with `Seurat`, `SingleCellExperiment`, and
  `SpatialExperiment` objects.
- Core function
  [`RankMap()`](https://github.com/jinming-cheng/RankMap/reference/RankMap.md)
  provides a streamlined pipeline for preprocessing, model training, and
  prediction.
- Customizable preprocessing: top-K gene masking, optional binning,
  expression weighting, and scaling.
- Additional functions:
  - [`computeRankedMatrix()`](https://github.com/jinming-cheng/RankMap/reference/ComputeRankedMatrix.md)
    – generate ranked matrices
  - [`trainRankModel()`](https://github.com/jinming-cheng/RankMap/reference/TrainRankModel.md)
    – train multinomial GLM
  - [`predictRankModel()`](https://github.com/jinming-cheng/RankMap/reference/PredictRankModel.md)
    – apply trained model to query data
  - [`evaluatePredictionPerformance()`](https://github.com/jinming-cheng/RankMap/reference/EvaluatePredictionPerformance.md)
    – assess accuracy
- Optimized for large datasets with significantly faster runtime than
  `SingleR`, `Azimuth`, and `RCTD`.
