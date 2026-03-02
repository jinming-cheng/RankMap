# Sample Cells by Cell Type with Minimum Representation

This function samples a specified total number of cells from a metadata
table, ensuring that each cell type is represented by at least a minimum
number of cells.

## Usage

``` r
SampleCellsByType(cell_metadata, n_total_cells = 5000, min_per_type = 50)
```

## Arguments

- cell_metadata:

  A data frame containing at least two columns: `cell_id` and
  `cell_type`.

- n_total_cells:

  Total number of cells to sample. Default is `5000`.

- min_per_type:

  Minimum number of cells to sample from each `cell_type`. Default is
  `50`.

## Value

A data frame of sampled cells with at least `min_per_type` cells per
`cell_type`.

## Examples

``` r
meta <- data.frame(
    cell_id = paste0("cell", 1:10000),
    cell_type = sample(c("T", "B", "Mac"), 10000, replace = TRUE)
)
set.seed(42)
sampled <- SampleCellsByType(meta,
    n_total_cells = 1000,
    min_per_type = 50
)
```
