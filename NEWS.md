
# RankMap 0.99.0

*Initial submission to Bioconductor*

## New Features

- Fast, robust, and scalable reference-based cell type annotation using multinomial regression on sparse ranked expression matrices.
- Supports both **single-cell** and **spatial transcriptomics** data.
- Compatible with `Seurat`, `SingleCellExperiment`, and `SpatialExperiment` objects.
- Core function `RankMap()` provides a streamlined pipeline for preprocessing, model training, and prediction.
- Customizable preprocessing: top-K gene masking, optional binning, expression weighting, and scaling.
- Additional functions:
  - `ComputeRankedMatrix()` – generate ranked matrices
  - `TrainRankModel()` – train multinomial GLM
  - `PredictRankModel()` – apply trained model to query data
  - `EvaluatePredictionPerformance()` – assess accuracy
- Optimized for large datasets with significantly faster runtime than `SingleR`, `Azimuth`, and `RCTD`.

