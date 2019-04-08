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

realcounts <- files %>% select(name, url)
usethis::use_data(realcounts, compress = "xz", overwrite = TRUE)

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
