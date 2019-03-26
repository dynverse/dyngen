library(tidyverse)

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

output_dir <- tempfile()
dir.create(output_dir)

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

regulatorycircuits_01_neurons_fetal_brain <- read_rds(paste0(output_dir, "/regulatorycircuits_01_neurons_fetal_brain.rds"))

usethis::use_data(regulatorycircuits_01_neurons_fetal_brain, compress = "xz", overwrite = TRUE)

realnets <- metadf %>% select(name, url)
usethis::use_data(realnets, compress = "xz", overwrite = TRUE)
