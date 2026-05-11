
# RankMap 1.1.1

- Add `topCellTypeGenes()` function

# RankMap 1.0.0

Bioconductor 3.23 Release

# RankMap 0.99.1

## Changes

- Updated `DESCRIPTION` to require `R (>= 4.6.0)`.
- Refactored repeated matrix transposition logic in `trainRankModel()` and `predictRankModel()` into an internal helper function (`.transposeMatrix`).
- Reduced use of `::` by importing selected functions (e.g., `glmnet()`, `cv.glmnet()`) to improve dependency checking.
- Added package-level documentation (`RankMap-package.R`).
- Updated vignette:
  - Removed GitHub installation instructions.
  - Added Bioconductor installation instructions using `BiocManager::install("RankMap")`.
  - Added unique labels to all code chunks.


# RankMap 0.99.0

*Initial submission to Bioconductor*

## New Features

- Fast, robust, and scalable reference-based cell type annotation using multinomial regression on ranked expression matrices.
- Supports both **single-cell** and **spatial transcriptomics** data.
- Compatible with `Seurat`, `SingleCellExperiment`, and `SpatialExperiment` objects.
- Core function `RankMap()` provides a streamlined pipeline for preprocessing, model training, and prediction.
- Customizable preprocessing: top-K gene masking, optional binning, expression weighting, and scaling.
- Additional functions:
  - `computeRankedMatrix()` – generate ranked matrices
  - `trainRankModel()` – train multinomial GLM
  - `predictRankModel()` – apply trained model to query data
  - `evaluatePredictionPerformance()` – assess accuracy
- Optimized for large datasets with significantly faster runtime than `SingleR`, `Azimuth`, and `RCTD`.

