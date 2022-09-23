#' Convert simulation output to different formats.
#' 
#' For use with other packages compatible with dyno, anndata, SingleCellExperiment, or Seurat.
#' 
#' @param model A dyngen output model for which the experiment has been emulated with [generate_experiment()].
#' @param store_cellwise_grn Whether or not to also store cellwise GRN information.
#' @param store_dimred Whether or not to store the dimensionality reduction constructed on the true counts.
#' @param store_rna_velocity WHether or not to store the log propensity ratios.
#' 
#' @return A dataset object.
#' 
#' @export
#' @rdname convert
#' 
#' @examples
#' data("example_model")
#' dataset <- wrap_dataset(example_model)
as_dyno <- function(
  model,
  store_dimred = !is.null(model$simulations$dimred),
  store_cellwise_grn = !is.null(model$experiment$cellwise_grn),
  store_rna_velocity = !is.null(model$experiment$rna_velocity)
) {
  requireNamespace("dynwrap")
  
  assert_that(
    !is.null(model$experiment), 
    msg = "model should be an object that was initialised with `initialise_model()`."
  )
  
  counts <- model$experiment$counts_mrna + model$experiment$counts_premrna
  dataset <- dynwrap::wrap_expression(
    id = model$id,
    counts = counts,
    counts_spliced = model$experiment$counts_mrna,
    counts_unspliced = model$experiment$counts_premrna,
    counts_protein = model$experiment$counts_protein,
    expression = as(log2(counts + 1), "CsparseMatrix"),
    cell_info = model$experiment$cell_info %>% select(-.data$from, -.data$to, -.data$time),
    feature_info = model$experiment$feature_info
  ) %>% 
    dynwrap::add_trajectory(
      milestone_network = model$gold_standard$network,
      progressions = model$experiment$cell_info %>%
        select(.data$cell_id, .data$from, .data$to, percentage = .data$time)
    )
  
  if (store_dimred) {
    dimred <- model$simulations$dimred[model$experiment$cell_info$step_ix, , drop = FALSE]
    rownames(dimred) <- model$experiment$cell_info$cell_id
    
    dataset <- dataset %>% 
      dynwrap::add_dimred(
        dimred = dimred,
        dimred_segment_points = model$gold_standard$dimred[!model$gold_standard$meta$burn, , drop = FALSE],
        dimred_segment_progressions = model$gold_standard$meta[!model$gold_standard$meta$burn, , drop = FALSE] %>% 
          select(.data$from, .data$to, percentage = .data$time)
      )
  }
  
  # add a few more values
  if (store_cellwise_grn) {
    regulatory_network <- model$feature_network %>% 
      select(regulator = .data$from, target = .data$to, .data$strength, .data$effect)
    regulation_sc <- model$experiment$cellwise_grn
    
    regulators <- unique(regulatory_network$regulator)
    targets <- colnames(dataset$counts)
    regulatory_network_sc <- 
      Matrix::summary(regulation_sc) %>% 
      transmute(
        cell_id = factor(dataset$cell_ids[.data$i], levels = dataset$cell_ids),
        regulator = factor(regulatory_network$regulator[.data$j], levels = regulators),
        target = factor(regulatory_network$target[.data$j], levels = targets),
        strength = .data$x
      ) %>% 
      as_tibble()
    
    
    dataset <- dataset %>% 
      dynwrap::add_regulatory_network(
        regulatory_network = regulatory_network,
        regulatory_network_sc = regulatory_network_sc,
        regulators = regulators,
        targets = targets
      )
  }
  
  if (store_rna_velocity) {
    dataset <- dataset %>% extend_with(
      "dynwrap::with_rna_velocity",
      rna_velocity = model$experiment$rna_velocity
    )
  }
  
  dataset
}


#' @rdname convert
#' @importFrom tibble column_to_rownames
#' @export
as_anndata <- function(
  model,
  store_dimred = !is.null(model$simulations$dimred),
  store_cellwise_grn = !is.null(model$experiment$cellwise_grn),
  store_rna_velocity = !is.null(model$experiment$rna_velocity)
) {
  requireNamespace("anndata")
  
  assert_that(
    !is.null(model$experiment), 
    msg = "model should be an object that was initialised with `initialise_model()`."
  )
  
  counts <- model$experiment$counts_mrna + model$experiment$counts_premrna
  
  ad_args <- list(
    X = counts,
    layers = list(
      counts_spliced = model$experiment$counts_mrna,
      counts_unspliced = model$experiment$counts_premrna,
      counts_protein = model$experiment$counts_protein,
      logcounts = as(log2(counts + 1), "CsparseMatrix")
    ),
    obs = model$experiment$cell_info %>%
      select(-.data$from, -.data$to, -.data$time) %>% 
      as.data.frame() %>% 
      column_to_rownames("cell_id"),
    var = model$experiment$feature_info %>% 
      as.data.frame() %>% 
      column_to_rownames("feature_id"),
    uns = list(
      traj_milestone_network = model$gold_standard$network,
      traj_progressions = model$experiment$cell_info %>%
        select(.data$cell_id, .data$from, .data$to, percentage = .data$time) %>% 
        as.data.frame()
    )
  )
  
  if (store_dimred) {
    dimred <- model$simulations$dimred[model$experiment$cell_info$step_ix, , drop = FALSE]
    rownames(dimred) <- model$experiment$cell_info$cell_id
    ad_args$obsm$dimred <- dimred
    
    ad_args$uns$traj_dimred_segments <- bind_cols(
      model$gold_standard$meta[!model$gold_standard$meta$burn, , drop = FALSE],
      as.data.frame(model$gold_standard$dimred[!model$gold_standard$meta$burn, , drop = FALSE])
    )
  }
  
  if (store_cellwise_grn) {
    ad_args$uns$regulatory_network <- model$feature_network %>% 
      select(regulator = .data$from, target = .data$to, .data$strength, .data$effect)
    ad_args$obsm$regulatory_network_sc <- model$experiment$cellwise_grn
    ad_args$uns$regulatory_network_regulators <- unique(model$feature_network$from)
    ad_args$uns$regulatory_network_targets <- colnames(counts)
  }
  
  if (store_rna_velocity) {
    ad_args$layers$rna_velocity <- model$experiment$rna_velocity
  }
  
  do.call(anndata::AnnData, ad_args)
}

#' @rdname convert
#' @importFrom tibble column_to_rownames
#' @export
as_sce <- function(
  model,
  store_dimred = !is.null(model$simulations$dimred),
  store_cellwise_grn = !is.null(model$experiment$cellwise_grn),
  store_rna_velocity = !is.null(model$experiment$rna_velocity)
) {
  requireNamespace("SingleCellExperiment")
  requireNamespace("SummarizedExperiment")
  
  assert_that(
    !is.null(model$experiment), 
    msg = "model should be an object that was initialised with `initialise_model()`."
  )
  
  counts <- model$experiment$counts_mrna + model$experiment$counts_premrna
  
  sce_args <- list(
    assays = list(
      counts = t(counts),
      logcounts = t(as(log2(counts + 1), "CsparseMatrix")),
      counts_spliced = t(model$experiment$counts_mrna),
      counts_unspliced = t(model$experiment$counts_premrna),
      counts_protein = t(model$experiment$counts_protein)
    ),
    colData = model$experiment$cell_info %>%
      select(-.data$from, -.data$to, -.data$time) %>% 
      as.data.frame() %>% 
      column_to_rownames("cell_id"),
    rowData = model$experiment$feature_info %>% 
      as.data.frame() %>% 
      column_to_rownames("feature_id"),
    metadata = list(
      traj_milestone_network = model$gold_standard$network,
      traj_progressions = model$experiment$cell_info %>%
        select(.data$cell_id, .data$from, .data$to, percentage = .data$time) %>% 
        as.data.frame()
    ),
    reducedDims = list(),
    altExps = list()
  )
  
  if (store_dimred) {
    dimred <- model$simulations$dimred[model$experiment$cell_info$step_ix, , drop = FALSE]
    rownames(dimred) <- model$experiment$cell_info$cell_id
    sce_args$reducedDims$MDS <- dimred
    
    sce_args$metadata$traj_dimred_segments <- bind_cols(
      model$gold_standard$meta[!model$gold_standard$meta$burn, , drop = FALSE],
      as.data.frame(model$gold_standard$dimred[!model$gold_standard$meta$burn, , drop = FALSE])
    )
  }
  
  if (store_cellwise_grn) {
    sce_args$altExps$regulatory_network <- SummarizedExperiment::SummarizedExperiment(
      assays = t(model$experiment$cellwise_grn),
      rowData = model$feature_network %>% 
        select(regulator = .data$from, target = .data$to, .data$strength, .data$effect),
      metadata = list(
        regulators = unique(model$feature_network$from),
        targets = colnames(counts)
      )
    )
  }
  
  if (store_rna_velocity) {
    sce_args$assays$rna_velocity <- t(model$experiment$rna_velocity)
  }
  
  do.call(SingleCellExperiment::SingleCellExperiment, sce_args)
}

#' @rdname convert
#' @importFrom tibble column_to_rownames
#' @export
as_seurat <- function(
  model,
  store_dimred = !is.null(model$simulations$dimred),
  store_cellwise_grn = !is.null(model$experiment$cellwise_grn),
  store_rna_velocity = !is.null(model$experiment$rna_velocity)
) {
  requireNamespace("Seurat")
  
  assert_that(
    !is.null(model$experiment), 
    msg = "model should be an object that was initialised with `initialise_model()`."
  )
  
  # collect metadata  
  feat_metadata <-
    model$experiment$feature_info %>% 
    mutate(feature_id = gsub("_", "-", .data$feature_id)) %>% # rename because seurat can't deal with underscores
    as.data.frame() %>% 
    column_to_rownames("feature_id")
  
  cell_metadata <- 
    model$experiment$cell_info %>%
    select(-.data$from, -.data$to, -.data$time) %>%
    as.data.frame() %>%
    column_to_rownames("cell_id")

  # create assays
  counts <- t(model$experiment$counts_mrna + model$experiment$counts_premrna)
  rownames(counts) <- rownames(feat_metadata)
  assay_obj <- Seurat::CreateAssayObject(counts = counts) %>% 
    Seurat::AddMetaData(feat_metadata)
  # logcounts = t(as(log2(counts + 1), "CsparseMatrix")),
  
  counts_mrna <- t(model$experiment$counts_mrna)
  rownames(counts_mrna) <- rownames(feat_metadata)
  counts_premrna <- t(model$experiment$counts_premrna)
  rownames(counts_premrna) <- rownames(feat_metadata)
  counts_protein <- t(model$experiment$counts_protein)
  rownames(counts_protein) <- rownames(feat_metadata)
  
  # construct seurat object and add assays
  seurat_obj <- suppressWarnings(Seurat::CreateSeuratObject(assay_obj, meta.data = cell_metadata))
  seurat_obj[["spliced"]] <- Seurat::CreateAssayObject(counts_mrna)
  seurat_obj[["unspliced"]] <- Seurat::CreateAssayObject(counts_premrna)
  seurat_obj[["protein"]] <- Seurat::CreateAssayObject(counts_protein)
  
  # add trajectory info
  Seurat::Misc(seurat_obj, "traj_milestone_network") <- model$gold_standard$network
  Seurat::Misc(seurat_obj, "traj_progressions") <- model$experiment$cell_info %>%
    select(.data$cell_id, .data$from, .data$to, percentage = .data$time) %>% 
    as.data.frame()
  
  if (store_dimred) {
    dimred <- model$simulations$dimred[model$experiment$cell_info$step_ix, , drop = FALSE]
    rownames(dimred) <- model$experiment$cell_info$cell_id
    dr <- Seurat::CreateDimReducObject(
      embeddings = dimred,
      assay = "RNA"
    )
    
    Seurat::Misc(dr, "traj_segments") <- bind_cols(
      model$gold_standard$meta[!model$gold_standard$meta$burn, , drop = FALSE],
      as.data.frame(model$gold_standard$dimred[!model$gold_standard$meta$burn, , drop = FALSE])
    )
    
    seurat_obj@reductions[["MDS"]] <- dr
  }
  
  if (store_cellwise_grn) {
    grn_feat_metadata <-
      model$feature_network %>% 
      select(regulator = .data$from, target = .data$to, .data$strength, .data$effect) %>% 
      mutate(row_name = gsub("_", "-", paste0(.data$regulator, "->", .data$target))) %>% 
      column_to_rownames("row_name")
    
    grn_cellwise <- t(model$experiment$cellwise_grn)
    rownames(grn_cellwise) <- rownames(grn_feat_metadata)
    
    grn_assay <- Seurat::CreateAssayObject(counts = grn_cellwise) %>% 
      Seurat::AddMetaData(grn_feat_metadata)
    Seurat::Misc(grn_assay, "regulators") <- unique(model$feature_network$from)
    Seurat::Misc(grn_assay, "targets") <- rownames(counts)
    
    seurat_obj[["regulatorynetwork"]] <- grn_assay
  }
  
  if (store_rna_velocity) {
    rna_velocity <- t(model$experiment$rna_velocity)
    rownames(rna_velocity) <- rownames(feat_metadata)
    seurat_obj[["rnavelocity"]] <- Seurat::CreateAssayObject(rna_velocity)
  }
  
  seurat_obj
}

convert_progressions_to_milestone_percentages <- function (progressions) {
  selfs <- progressions %>% filter(.data$from == .data$to) %>% select(.data$cell_id, milestone_id = .data$from) %>% mutate(percentage = 1)
  progressions <- progressions %>% filter(.data$from != .data$to)
  from_mls <- tapply(progressions$from, progressions$cell_id, first, default = NA_character_)
  from_pct <- 1 - tapply(progressions$percentage, progressions$cell_id, sum, default = NA_real_)
  froms <- tibble(
    cell_id = names(from_mls) %||% character(), 
    milestone_id = from_mls[.data$cell_id] %>% unname() %>% as.character(), 
    percentage = from_pct[.data$cell_id] %>% unname() %>% as.numeric()
  )
  tos <- progressions %>% select(.data$cell_id, milestone_id = .data$to, .data$percentage)
  bind_rows(selfs, froms, tos)
}

#' @rdname convert
#' @export
as_list <- function(
  model,
  store_dimred = !is.null(model$simulations$dimred),
  store_cellwise_grn = !is.null(model$experiment$cellwise_grn),
  store_rna_velocity = !is.null(model$experiment$rna_velocity)
) {
  assert_that(
    !is.null(model$experiment), 
    msg = "model should be an object that was initialised with `initialise_model()`."
  )
  
  counts <- model$experiment$counts_mrna + model$experiment$counts_premrna
  progressions <- model$experiment$cell_info %>%
    select(.data$cell_id, .data$from, .data$to, percentage = .data$time)
  dataset <- list(
    id = model$id,
    cell_ids = rownames(counts),
    feature_ids = colnames(counts),
    counts = counts,
    counts_spliced = model$experiment$counts_mrna,
    counts_unspliced = model$experiment$counts_premrna,
    counts_protein = model$experiment$counts_protein,
    expression = as(log2(counts + 1), "CsparseMatrix"),
    cell_info = model$experiment$cell_info %>% select(-.data$from, -.data$to, -.data$time),
    feature_info = model$experiment$feature_info,
    milestone_ids = unique(c(model$gold_standard$network$from, model$gold_standard$network$to)),
    milestone_network = model$gold_standard$network,
    progressions = progressions,
    milestone_percentages = convert_progressions_to_milestone_percentages(progressions)
  )
  
  if (store_dimred) {
    dimred <- model$simulations$dimred[model$experiment$cell_info$step_ix, , drop = FALSE]
    rownames(dimred) <- model$experiment$cell_info$cell_id
    
    dataset$dimred <- dimred
    dataset$dimred_segment_points <- model$gold_standard$dimred[!model$gold_standard$meta$burn, , drop = FALSE]
    dataset$dimred_segment_progressions <- model$gold_standard$meta[!model$gold_standard$meta$burn, , drop = FALSE] %>% 
      select(.data$from, .data$to, percentage = .data$time)
  }
  
  # add a few more values
  if (store_cellwise_grn) {
    regulatory_network <- model$feature_network %>% 
      select(regulator = .data$from, target = .data$to, .data$strength, .data$effect)
    regulation_sc <- model$experiment$cellwise_grn
    
    regulators <- unique(regulatory_network$regulator)
    targets <- colnames(dataset$counts)
    regulatory_network_sc <- 
      Matrix::summary(regulation_sc) %>% 
      transmute(
        cell_id = factor(dataset$cell_ids[.data$i], levels = dataset$cell_ids),
        regulator = factor(regulatory_network$regulator[.data$j], levels = regulators),
        target = factor(regulatory_network$target[.data$j], levels = targets),
        strength = .data$x
      ) %>% 
      as_tibble()
    
    
    dataset$regulatory_network <- regulatory_network
    dataset$regulatory_network_sc <- regulatory_network_sc
    dataset$regulators <- regulators
    dataset$targets <- targets
  }
  
  if (store_rna_velocity) {
    dataset$rna_velocity <- model$experiment$rna_velocity
  }
  
  dataset
}

conversion_funs <- list(
  dyno = as_dyno,
  sce = as_sce,
  seurat = as_seurat,
  anndata = as_anndata,
  list = as_list,
  none = function(...) NULL
)

#' @rdname convert
#' @param format Which output format to use, must be one of 'dyno' (requires `dynwrap`), 'sce' (requires `SingleCellExperiment`), 'seurat' (requires `Seurat`), 'anndata' (requires `anndata`), 'list' or 'none'.
#' @export
wrap_dataset <- function(
  model,
  format = c("list", "dyno", "sce", "seurat", "anndata", "none"),
  store_dimred = !is.null(model$simulations$dimred),
  store_cellwise_grn = !is.null(model$experiment$cellwise_grn),
  store_rna_velocity = !is.null(model$experiment$rna_velocity)
) {
  format <- match.arg(format)
  fun <- conversion_funs[[format]]
  
  fun(
    model = model, 
    store_dimred = store_dimred,
    store_cellwise_grn = store_cellwise_grn,
    store_rna_velocity = store_rna_velocity
  )
}