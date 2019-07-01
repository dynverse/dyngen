#' @importFrom utils download.file
.download_cacheable_file <- function(url, cache_dir, verbose) {
  file <- 
    if (is.null(cache_dir)) {
      fil <- tempfile()
      on.exit(file.remove(fil))
      fil
    } else {
      file <- paste0(
        sub("/*$", "/", cache_dir),
        sub(".*/", "", url)
      )
    }
  
  if (!file.exists(file)) {
    if (!dir.exists(dirname(file))) {
      dir.create(dirname(file), recursive = TRUE)
    }
    download.file(url, destfile = file, quiet = !verbose)
  }
  
  readRDS(file)
}