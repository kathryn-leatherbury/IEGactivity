# BEGIN: make expression-score scaling configurable while retaining the default factor of 100.
#' Compute per-nucleus immediate early gene scores
#'
#' Calculates IEG scores from expression or count data in a Seurat object.
#'
#' @param object Seurat object.
#' @param iegs Character vector of IEG identifiers to score.
#' @param gex.assay Assay to use. Defaults to "RNA".
#' @param cluster.level A factor in object metadata defining cell-type or cluster
#'   levels for variance estimation.
#' @param sample.id A factor or character in object metadata containing sample
#'   identifiers, used to give samples equal weight in variance estimation.
#' @param group.by A factor or character in object metadata to group nuclei within
#'   each cell type for variance estimation. NULL uses the cell type alone.
#' @param scoring.method Either "expression" or "count". Expression uses normalized gene expression scores. Count uses a sum of IEG counts per nucleus.
#' @param norm.weights Default TRUE. Apply weighted normalization using se variance-based gene weights for
#'   expression scoring; if FALSE, weight genes equally.
#' @param scale.factor Numeric multiplier applied to expression scores. Defaults
#'   to 100. Not used for count scoring.
#'
#' @return A named numeric vector with one score per nucleus in the selected
#'   assay layer.
#'
#' @examples
#' \dontrun{
#' seurat_obj$IEGscore <- ComputeIEGscore(
#'   object = seurat_obj,
#'   iegs = ieg_genes,
#'   cluster.level = "cell_type",
#'   sample.id = "sample_id"
#' )
#' }
#'
#' @export
ComputeIEGscore <- function(
    object,
    iegs,
    gex.assay      = 'RNA',
    cluster.level,
    sample.id,
    group.by       = NULL,
    scoring.method = 'expression',  # options either expression or count
    norm.weights   = TRUE, # normalize across samples by applying equal weights and standardize by group level if group.by is set
    scale.factor   = 100
) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject is required.")
  if (!scoring.method %in% c("count", "expression")) {
    stop("scoring.method must be either 'count' or 'expression'.")
  }
  if (!is.character(iegs) || !length(iegs) || anyNA(iegs) ||
      anyDuplicated(iegs) || any(!nzchar(iegs))) {
    stop("iegs must contain unique, nonempty gene names.")
  }
  md <- object@meta.data
  fields <- c(sample.id, group.by)
  if (scoring.method == "expression" && norm.weights == TRUE) {
    fields <- c(cluster.level, fields)
  }
  if (!all(fields %in% names(md))) {
    stop("Required metadata columns are missing: ",
         paste(setdiff(fields, names(md)), collapse = ", "))
  }
  samples <- unique(md[c(sample.id, group.by)])
  if (anyDuplicated(samples[[sample.id]])) stop("Each sample.id must map to only one group.by value.")

  if (scoring.method == "expression" && norm.weights == TRUE) {
    score.layer <- 'data'
    if (!score.layer %in% SeuratObject::Layers(object, assay = gex.assay, search = NA)) {
      stop("The exact requested layer is absent; resolve split layers first.")
    }
    available <- rownames(SeuratObject::LayerData(object, assay = gex.assay, layer = score.layer))
    missing <- setdiff(iegs, available)
    if (length(missing)) stop("Missing requested genes: ", paste(missing, collapse = ", "))

    mat <- as.matrix(SeuratObject::LayerData(object, assay = gex.assay, layer = score.layer)[
      iegs, , drop = FALSE])
    if (!ncol(mat) || any(!is.finite(mat)) || any(mat < 0)) {
      stop("Use a nonempty matrix of finite, nonnegative expression values.")
    }
    if (!all(colnames(mat) %in% rownames(md))) {
      stop("Required metadata or matching cell names are missing.")
    }
    md <- md[colnames(mat), fields, drop = FALSE]
    md[] <- lapply(md, as.character)

    # Calculate weights within each cell type when no grouping factor is supplied.
    if (is.null(group.by)) {
      names(md) <- c("celltype", "individual")
      if (anyNA(md) || any(vapply(md, function(z) any(!nzchar(trimws(z))), logical(1)))) {
        stop("Required metadata contain missing or blank entries.")
      }
      G <- nrow(mat)
      score <- stats::setNames(numeric(ncol(mat)), colnames(mat))
      for (ct in unique(md$celltype)) {
        ct_cells <- which(md$celltype == ct)
        ids <- unique(md$individual[ct_cells])
        means <- within_variance <- matrix(0, G, length(ids))
        for (i in seq_along(ids)) {
          j <- ct_cells[md$individual[ct_cells] == ids[i]]
          x <- mat[, j, drop = FALSE]
          means[, i] <- rowMeans(x)
          within_variance[, i] <- rowMeans(sweep(x, 1, means[, i], "-")^2)
        }
        # Equal sample contributions: within-sample plus between-sample variation.
        celltype_mean <- rowMeans(means)
        variance <- rowMeans(within_variance) +
          rowMeans(sweep(means, 1, celltype_mean, "-")^2)
        positive <- variance[variance > 0]
        stabilizer <- if (length(positive)) stats::median(positive) else 0
        if (stabilizer > 0) {
          adjusted_variance <- variance + stabilizer
          w <- 1 / sqrt(adjusted_variance)
          w <- w / sum(w)
        } else {
          w <- rep(1 / G, G)
        }
        score[ct_cells] <- as.numeric(w %*% mat[, ct_cells, drop = FALSE]) * scale.factor
      }
      return(score)
    }

    names(md) <- c("celltype", "individual", "condition")
    if (anyNA(md) || any(vapply(md, function(z) any(!nzchar(trimws(z))), logical(1)))) {
      stop("Required metadata contain missing or blank entries.")
    }
    # RUN IF FUNCTION IS SET TO EXPRESSION AND norm.weights IS SET TO TRUE
    G <- nrow(mat)
    score <- stats::setNames(numeric(ncol(mat)), colnames(mat))
    for (ct in unique(md$celltype)) {
      ct_cells <- which(md$celltype == ct)
      states <- unique(md$condition[ct_cells])
      status_cells <- lapply(states, function(s) ct_cells[md$condition[ct_cells] == s])
      variance <- matrix(0, G, length(states))

      for (s in seq_along(states)) {
        cells <- status_cells[[s]]
        ids <- unique(md$individual[cells])
        means <- within_variance <- matrix(0, G, length(ids))
        for (i in seq_along(ids)) {
          j <- cells[md$individual[cells] == ids[i]]
          x <- mat[, j, drop = FALSE]
          means[, i] <- rowMeans(x)
          within_variance[, i] <- rowMeans(sweep(x, 1, means[, i], "-")^2)
        }
        # Equal sample contributions: within-sample plus between-sample variation.
        status_mean <- rowMeans(means)
        variance[, s] <- rowMeans(within_variance) +
          rowMeans(sweep(means, 1, status_mean, "-")^2)
      }

      # Equal group contributions; differences between group means are excluded.
      reference <- rowMeans(variance)
      positive <- reference[reference > 0]
      stabilizer <- if (length(positive)) stats::median(positive) else 0

      for (s in seq_along(states)) {
        if (stabilizer > 0) {
          adjusted_variance <- 0.5 * variance[, s] + 0.5 * reference + stabilizer
          w <- 1 / sqrt(adjusted_variance)
          w <- w / sum(w)
        } else {
          w <- rep(1 / G, G)
        }
        j <- status_cells[[s]]
        score[j] <- as.numeric(w %*% mat[, j, drop = FALSE]) * scale.factor
      }
    }
  } else if (scoring.method == "expression" && norm.weights == FALSE) {
    # RUN IF FUNCTION IS SET TO EXPRESSION AND norm.weights IS SET TO FALSE
    score.layer <- 'data'
    if (!score.layer %in% SeuratObject::Layers(object, assay = gex.assay, search = NA)) {
      stop("The exact requested layer is absent; resolve split layers first.")
    }
    available <- rownames(SeuratObject::LayerData(object, assay = gex.assay, layer = score.layer))
    missing <- setdiff(iegs, available)
    if (length(missing)) stop("Missing requested genes: ", paste(missing, collapse = ", "))
    invisible(utils::capture.output(
      ## read this score input from gex.assay.
      mat <- as.matrix(SeuratObject::LayerData(object, assay = gex.assay, layer = score.layer)[
        iegs, , drop = FALSE])
      ## retain the package argument name assay.
    ))
    n_ieg <- length(iegs)
    score <- (Matrix::colSums(mat, na.rm = TRUE)/n_ieg) * scale.factor
    score[score == 'NaN'] = 0
  } else if (scoring.method == 'count'){
    # RUN IF FUNCTION IS SET TO COUNT
    score.layer <- 'counts'
    if (!score.layer %in% SeuratObject::Layers(object, assay = gex.assay, search = NA)) {
      stop("The exact requested layer is absent; resolve split layers first.")
    }
    available <- rownames(SeuratObject::LayerData(object, assay = gex.assay, layer = score.layer))
    missing <- setdiff(iegs, available)
    if (length(missing)) stop("Missing requested genes: ", paste(missing, collapse = ", "))
    invisible(utils::capture.output(
      ## read this score input from gex.assay.
      mat <- as.matrix(SeuratObject::LayerData(object, assay = gex.assay, layer = score.layer)[
        iegs, , drop = FALSE])
      ## retain the package argument name assay.
    ))
    mat[mat > 0] <- 1
    score <- Matrix::colSums(mat)
  }
  score
}
# END: configurable expression-score scaling.
