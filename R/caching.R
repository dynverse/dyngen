.download_cacheable_file <- function(url, model) {
  file <- 
    if (is.null(model$download_cache_dir)) {
      fil <- tempfile()
      on.exit(file.remove(fil))
      fil
    } else {
      file <- paste0(
        sub("/*$", "/", model$download_cache_dir),
        sub(".*/", "", url)
      )
    }
  
  if (!file.exists(file)) {
    if (!dir.exists(dirname(file))) {
      dir.create(dirname(file), recursive = TRUE)
    }
    download.file(url, destfile = file, quiet = !model$verbose)
  }
  
  readRDS(file)
}