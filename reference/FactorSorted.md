# Order factor levels by frequency

Takes a vector, converts it to a factor, and reorders its levels by
their frequency of occurrence (most to least frequent by default).

## Usage

``` r
factorSorted(x, decreasing = TRUE)
```

## Arguments

- x:

  A vector (character or factor) to be reordered by level frequency.

- decreasing:

  Logical. If TRUE (default), levels are ordered from most to least
  frequent; if FALSE, from least to most frequent.

## Value

A factor with levels ordered by frequency.

## Details

Frequencies are computed with `table(x)`, which ignores `NA` values by
default. `NA` entries in `x` are preserved in the output but are not
considered a level.

## Examples

``` r
factorSorted(c("a", "b", "a", "c", "b", "a"))
#> [1] a b a c b a
#> Levels: a b c
factorSorted(c("a", "b", "a", "c", "b", "a"), decreasing = FALSE)
#> [1] a b a c b a
#> Levels: c b a
factorSorted(factor(c("x", "y", "x", NA)))
#> [1] x    y    x    <NA>
#> Levels: x y
```
