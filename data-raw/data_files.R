library(tidyverse)


output_dir <- tempfile()
dir.create(output_dir)

#################################################
##              DOWNLOAD REALNETS              ##
#################################################
tmp_zip <- tempfile(fileext = ".zip")
download.file("http://www2.unil.ch/cbg/regulatorycircuits/Network_compendium.zip", tmp_zip)
on.exit(file.remove(tmp_zip))

tmp_dir <- tempfile()
unzip(tmp_zip, exdir = tmp_dir)
on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE))

metadf <- 
  list.files(
    paste0(tmp_dir, "/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks"), 
    full.names = TRUE,
    pattern = "txt.gz"
  ) %>% 
  tibble(file = .) %>% 
  mutate(
    name = file %>% str_replace("^.*/", "") %>% str_replace("\\.txt\\.gz$", "") %>% paste0("regulatorycircuits_", .),
    url = paste0("https://github.com/dynverse/dyngen/raw/data_files/", name, ".rds")
  )

pbapply::pblapply(
  seq_len(nrow(metadf)),
  cl = 8,
  function(i) {
    file <- metadf$file[[i]]
    name <- metadf$name[[i]]
    
    df <- 
      read_tsv(file, col_names = c("regulator", "target", "weight"), col_types = c("regulator" = "c", "target" = "c", "weight" = "d")) %>% 
      arrange(desc(weight)) %>% 
      slice(seq_len(250000))
    
    genes <- unique(c(df$regulator, df$target))
    
    m <- Matrix::sparseMatrix(
      i = match(df$regulator, genes),
      j = match(df$target, genes),
      x = df$weight,
      dimnames = list(
        genes[sort(unique(match(df$regulator, genes)))],
        genes
      )
    )
    
    write_rds(m, paste0(output_dir, "/", name, ".rds"), compress = "xz")
  }
)

cat("Output is stored at ", output_dir, "\n", sep = "")

# regulatorycircuits_01_neurons_fetal_brain <- read_rds(paste0(output_dir, "/regulatorycircuits_01_neurons_fetal_brain.rds"))
# usethis::use_data(regulatorycircuits_01_neurons_fetal_brain, compress = "xz", overwrite = TRUE)

realnets <- metadf %>% select(name, url)
usethis::use_data(realnets, compress = "xz", overwrite = TRUE)

#################################################
##              DOWNLOAD DYNREAL               ##
#################################################

library(httr)

# config
deposit_id <- 1443566

# retrieve file metadata from zenodo
files <-
  GET(glue::glue("https://zenodo.org/api/records/{deposit_id}")) %>%
  httr::content() %>%
  .$files %>%
  map_df(function(l) {
    as_tibble(t(unlist(l)))
  }) %>%
  filter(grepl("^real/", filename)) %>% 
  mutate(
    name = filename %>% str_replace_all("\\.rds$", "") %>% str_replace_all("/", "_") %>% paste0("zenodo_", deposit_id, "_", .),
    url = paste0("https://github.com/dynverse/dyngen/raw/data_files/", name, ".rds"),
    local_out = paste0(output_dir, "/", name, ".rds")
  )

# iterate over rows
pbapply::pblapply(
  seq_len(nrow(files)),
  cl = 8,
  function(i) {
    tmp <- tempfile()
    on.exit(file.remove(tmp))
    
    # download file
    download.file(files$links.download[[i]], tmp, quiet = TRUE)
    
    # read counts
    counts <- read_rds(tmp)$counts %>% Matrix::Matrix(sparse = TRUE)
    
    # write data
    write_rds(counts, files$local_out[[i]], compress = "xz")
    
    invisible()
  }
)

## Add one more dataset
new_file <- tibble(
  name = "GSE100866_CBMC_8K_13AB_10X-RNA_umi",
  url = paste0("https://github.com/dynverse/dyngen/raw/data_files/", name, ".rds"),
  local_out = paste0(output_dir, "/", name, ".rds")
)

if (!file.exists(new_file$local_out)) {
  tmp <- tempfile()
  download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", tmp, quiet = TRUE)
  cbmc.rna <- read.csv(file = tmp, sep = ",", header = TRUE, row.names = 1)
  cbmc_counts <- Matrix::Matrix(as.matrix(cbmc.rna), sparse = TRUE) %>% Matrix::t()
  cells_expressed <- Matrix::colMeans(cbmc_counts > 0)
  cbmc_counts <- cbmc_counts[, cells_expressed > .05]
  
  write_rds(cbmc_counts, new_file$local_out, compress = "xz")
}
files <- bind_rows(files, new_file)



# perform checks and filtering
# files %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
realcounts <- pmap_df(files, function(name, url, local_out, ...) {
  cou <- readr::read_rds(local_out)
  
  library_size <- Matrix::rowSums(cou)
  
  cpm <- sweep(cou, 1, 1e6 / library_size, "*")
  
  num_cells <- nrow(cou)
  num_features <- ncol(cou)
  
  dropout_rate <- mean(as.matrix(cou) == 0)
  average_log2_cpm <- apply(cpm, 2, mean)
  variance_log2_cpm <- apply(cpm, 2, var)
  
  tibble(
    name, 
    url,
    num_cells,
    num_features,
    median_library_size = median(library_size),
    dropout_rate,
    median_average_log2_cpm = median(average_log2_cpm),
    median_variance_log2_cpm = median(variance_log2_cpm)
  )
})

ggplot(realcounts) + 
  geom_point(aes(num_cells, median_library_size, alpha = num_cells > 300 & median_library_size > 1000, colour = dropout_rate)) + 
  scale_y_log10() + scale_x_log10() + viridis::scale_color_viridis() + theme_bw()

qc_pass_names <- c(
  "zenodo_1443566_real_gold_cellbench-SC1_luyitian",
  "zenodo_1443566_real_gold_cellbench-SC2_luyitian",
  "zenodo_1443566_real_gold_cellbench-SC3_luyitian",
  "zenodo_1443566_real_gold_developing-dendritic-cells_schlitzer",
  "zenodo_1443566_real_gold_psc-astrocyte-maturation-neuron_sloan",
  "zenodo_1443566_real_silver_bone-marrow-mesenchyme-erythrocyte-differentiation_mca",
  "zenodo_1443566_real_silver_cell-cycle_leng",
  "zenodo_1443566_real_silver_embronic-mesenchyme-neuron-differentiation_mca",
  "zenodo_1443566_real_silver_embryonic-mesenchyme-stromal-cell-cxcl14-cxcl12-axis_mca",
  "zenodo_1443566_real_silver_mammary-gland-involution-endothelial-cell-aqp1-gradient_mca",
  "zenodo_1443566_real_silver_mouse-cell-atlas-combination-10",
  "zenodo_1443566_real_silver_mouse-cell-atlas-combination-4",
  "zenodo_1443566_real_silver_mouse-cell-atlas-combination-5",
  "zenodo_1443566_real_silver_mouse-cell-atlas-combination-8",
  "zenodo_1443566_real_silver_olfactory-projection-neurons-DA1_horns",
  "zenodo_1443566_real_silver_placenta-trophoblast-differentiation_mca",
  "zenodo_1443566_real_silver_thymus-t-cell-differentiation_mca",
  "zenodo_1443566_real_silver_trophoblast-stem-cell-trophoblast-differentiation_mca"
)
realcounts <- realcounts %>% mutate(qc_pass = name %in% qc_pass_names)

# ggplot(realcounts %>% mutate(qc_pass = row_number() %in% c(3, 4, 5, 6, 8, 24, 29, 30, 31, 34, 35, 51, 52, 56, 57, 60, 66, 72, 108, 110, 111))) + 
#   geom_point(aes(num_cells, median_library_size, alpha = qc_pass, colour = dropout_rate)) + 
#   scale_y_log10() + scale_x_log10() + viridis::scale_color_viridis() + theme_bw()
# 
# ggplot(realcounts %>% mutate(qc_pass = row_number() %in% c(3, 4, 5, 6, 8, 24, 29, 30, 31, 34, 35, 51, 52, 56, 57, 60, 66, 72, 108, 110, 111))) + 
#   geom_point(aes(median_average_log2_cpm, median_variance_log2_cpm, alpha = qc_pass, colour = dropout_rate)) + 
#   scale_y_log10() + scale_x_log10() + viridis::scale_color_viridis() + theme_bw()

usethis::use_data(realcounts, compress = "xz", overwrite = TRUE)


##################################################
##              DOWNLOAD KINETICS               ##
##################################################
# download data from Schwannhäusser et al., 2011, doi.org/10.1038/nature10098
tmpfile <- tempfile()
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fnature10098/MediaObjects/41586_2011_BFnature10098_MOESM304_ESM.xls", tmpfile)
tab <- readxl::read_excel(tmpfile) %>% 
  transmute(
    protein_ids = `Protein IDs`,
    protein_names = `Protein Names`,
    gene_names = `Gene Names`,
    uniprot_ids = `Uniprot IDs`,
    refset_protein_ids = `RefSeq protein IDs`,
    refseq_mrna_id = `Refseq mRNA ID`,
    ensembl_id = `ENSEMBL ID`,
    mgi_id = `MGI ID`,
    protein_length = as.numeric(`Protein length [amino acids]`),
    protein_copynumber = as.numeric(`Protein copy number average [molecules/cell]`),
    protein_halflife = as.numeric(`Protein half-life average [h]`),
    mrna_copynumber = as.numeric(`mRNA copy number average [molecules/cell]`),
    mrna_halflife = as.numeric(`mRNA half-life average [h]`),
    transcription_rate = as.numeric(`transcription rate (vsr) average [molecules/(cell*h)]`),
    translation_rate = as.numeric(`translation rate constant (ksp) average [molecules/(mRNA*h)]`)
  )

write_rds(tab, paste0(output_dir, "/schwannhausser2011.rds"), compress = "xz")

# get median values for deriving parameters for the TF backbone
tab %>% summarise_if(is.numeric, median, na.rm = TRUE)

# perform data imputation for sampling kinetics from this dataset
out <- mice::mice(tab %>% select_if(is.numeric), m = 1, maxit = 10, method = "pmm", seed = 1)

tab_imp <- mice::complete(out)

write_rds(tab_imp, paste0(output_dir, "/schwannhausser2011_imputed.rds"), compress = "xz")



#################################################
##    COMMIT DATA TO SEPARATE BRANCH    ##
#################################################
system(paste0(
  "pushd ", output_dir, "\n", 
  "git init\n",
  "git checkout --orphan data_files\n",
  "git add *.rds\n",
  "git commit --allow-empty -m 'build data files [ci skip]'\n",
  "git remote add origin git@github.com:dynverse/dyngen.git\n",
  "git push --force origin data_files\n",
  "popd\n"
))
