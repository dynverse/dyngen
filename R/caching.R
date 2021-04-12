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
    
    
    status <- suppressWarnings(
      utils::download.file(url, destfile = file, quiet = !verbose)
    )
    
    if (status != 0) stop("Cannot download file from ", url, call. = FALSE)
  }
  
  readRDS(file)
}