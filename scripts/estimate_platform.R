library(splatter)
library(Seurat)
library(tidyverse)

mutate <- dplyr::mutate
filter <- dplyr::filter
slice <- dplyr::slice

# dataset_id <- "psc-astrocyte-maturation-neuron_sloan"
dataset_id <- "hematopoiesis-clusters_olsson"
counts <- readRDS(glue::glue("../dynalysis/analysis/data/derived_data/datasets/real/{dataset_id}/original_counts.rds"))
dataset <- readRDS(glue::glue("../dynalysis/analysis/data/derived_data/datasets/real/{dataset_id}/dataset.rds"))
counts <- counts[dataset$cell_grouping$cell_id,]
counts %>% dim

# calculate changing genes
seurat <- Seurat::CreateSeuratObject(t(counts))
seurat@ident <- dataset$cell_grouping %>% slice(match(rownames(counts), cell_id)) %>% pull(group_id) %>% factor() %>% setNames(rownames(counts))
changing <- FindAllMarkers(seurat, logfc.treshold = 1, min.pct=0.4)
n_changing <- changing %>% filter(abs(avg_logFC) >= 1) %>% pull(gene) %>% unique() %>% length()
pct_changing <- n_changing / ncol(counts)

# data("sc_example_counts")
# counts <- t(sc_example_counts)

# splatter estimate libSizes etc
# splatter requires non-zero rows and columns to estimate
cell_sds <- counts %>% apply(1, sd)
gene_sds <- counts %>% apply(2, sd)
counts <- counts[cell_sds > 0, gene_sds > 0]
counts %>% dim

# estimate splatter parameters
estimate <- splatEstimate(t(counts))
class(estimate) <- "TheMuscularDogBlinkedQuietly." # change the class, so scater won't get magically loaded when the platform is loaded

platform <- lst(
  estimate,
  pct_changing,
  n_cells = nrow(counts),
  id = dataset_id
)

saveRDS(platform, glue::glue("inst/ext_data/platforms/{dataset_id}.rds"))
