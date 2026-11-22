# IEGactivity
Custom R functions for immediate early gene analyses in single cell multi-omics data. Compatible across multiple modalities for use with multiple conditions, sexes, batches and covariates. Package includes functions built for a standardized implementation of IEG analyses, including IEG modeling, IEG per nucleus scoring, and plotting visualizations

## Installation

Install the package from GitHub:

```r
install.packages("remotes")
remotes::install_github("kathryn-leatherbury/IEGactivity")
```

The package requires SeuratObject 5.0.0 or later and Matrix. Its other direct
dependencies, stats and utils, are included with R.

## ComputeIEGscore

`ComputeIEGscore()` returns a named vector of scores, with one value per nucleus
in the selected assay layer.

- **Weighted expression:** `scoring.method = "expression"` and
  `norm.weights = TRUE` (the defaults). Calculates a weighted mean of the
  requested IEG expression values, multiplied by 100.
- **Unweighted expression:** `scoring.method = "expression"` and
  `norm.weights = FALSE`. Calculates the mean of the requested IEG expression
  values, multiplied by 100.
- **Count:** `scoring.method = "count"`. Converts positive entries in the
  requested IEG count matrix to 1, then sums them for each nucleus.

Expression scoring reads the assay's `data` layer; count scoring reads its
`counts` layer. Gene identifiers must match the feature names in that layer.
The function checks for the exact layer name and reports missing genes.
Assays with split layers must have the requested layer resolved before scoring.

### Weighted scoring without grouping

Supply your Seurat object and a character vector of IEG identifiers. Replace
the example metadata names below with the names used in your object.

```r
seurat_obj$IEGscore <- IEGactivity::ComputeIEGscore(
  object = seurat_obj,
  iegs = ieg_genes,
  gex.assay = "RNA",
  cluster.level = "cell_type",
  sample.id = "sample_id"
)
```

`group.by` defaults to `NULL`. Gene variances and weights are calculated within
each cell type. Samples identified by `sample.id` contribute equally to variance
estimation.

### Weighted scoring with grouping

```r
seurat_obj$IEGscore <- IEGactivity::ComputeIEGscore(
  object = seurat_obj,
  iegs = ieg_genes,
  cluster.level = "cell_type",
  sample.id = "sample_id",
  group.by = "condition"
)
```

With `group.by` supplied, variances are calculated separately for its groups
within each cell type. The weighting calculation combines each group's gene
variance with the average gene variance across groups. Differences between
group means are excluded from that variance estimate. Each sample identifier
must map to one value of the supplied grouping variable.

### Unweighted expression and count scoring

```r
seurat_obj$IEGscore_unweighted <- IEGactivity::ComputeIEGscore(
  object = seurat_obj,
  iegs = ieg_genes,
  sample.id = "sample_id",
  scoring.method = "expression",
  norm.weights = FALSE
)

seurat_obj$IEGcount <- IEGactivity::ComputeIEGscore(
  object = seurat_obj,
  iegs = ieg_genes,
  sample.id = "sample_id",
  scoring.method = "count"
)
```

`cluster.level` is used for weighted expression scoring.
`norm.weights` does not affect count scoring.

For all parameter descriptions, see:

```r
help("ComputeIEGscore", package = "IEGactivity")
```

## License

MIT. See [LICENSE.md](LICENSE.md).
