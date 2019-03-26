library(tidyverse)

tmp_zip <- tempfile(fileext = ".zip")
download.file("http://www2.unil.ch/cbg/regulatorycircuits/Network_compendium.zip", tmp_zip)
on.exit(file.remove(tmp_zip))

tmp_dir <- tempfile()
unzip(tmp_zip, exdir = tmp_dir)
on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE))

set.seed(1)

network_files <- list.files(
  paste0(tmp_dir, "/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks"), 
  full.names = TRUE,
  pattern = "txt.gz"
) %>% 
  sample(4) %>% 
  sort()

network_names <- network_files %>% str_replace("^.*/", "") %>% str_replace("\\.txt\\.gz$", "")

realnets <- map(
  network_files,
  function(file) {
    df <- 
      read_tsv(file, col_names = c("regulator", "target", "weight"), col_types = c("regulator" = "c", "target" = "c", "weight" = "d")) %>% 
      arrange(desc(weight)) %>% 
      slice(seq_len(250000))
    
    genes <- unique(c(df$regulator, df$target))
    
    df %>% 
      mutate(
        regulator = match(regulator, genes),
        target = match(target, genes)
      )
  }
) %>% 
  set_names(network_names)

pryr::object_size(realnet)
usethis::use_data(realnet, compress = "xz", overwrite = TRUE)
