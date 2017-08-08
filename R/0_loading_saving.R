#' List all datasets
#' @import dplyr
#' @export
list_datasets = function() overviewer("datasets")

#' Save a dataset of a particular type
#' @import dplyr
#' @export
saver <- function(x, type) {
  newoverview <- map(x, ~.$info) %>% bind_rows()
  
  for(xi in x) {
    path <- paste0(.datasets_location, "/", type, "/", xi$info$id, ".rds")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(xi, path)
  }
  
  overview_location <- paste0(.datasets_location, "/", type, ".rds")
  overview <- if(file.exists(overview_location)) {readRDS(paste0(.datasets_location, "/", type, ".rds"))} else {tibble()  }
  
  overview %>% 
    filter(!(id %in% newoverview$id)) %>% 
    bind_rows(newoverview) %>% 
    saveRDS(paste0(.datasets_location, "/", type, ".rds"))
}

#' Generate an overview data.frame for a particular type of data
#' @import dplyr
#' @export
overviewer <- function(type) {
  overview_location <- paste0(.datasets_location, "/", type, ".rds")
  if(file.exists(overview_location)) {readRDS(paste0(.datasets_location, "/", type, ".rds"))} else {tibble(id=character())  }
}

#' Load multiple datasets of one type
#' @import dplyr
#' @export
loader <- function(x, type) {
  overview <- overviewer(type)
  
  map(x, function(xi) {
    path <- paste0(.datasets_location, "/", type, "/", xi, ".rds")
    readRDS(path)
  })
}






#' Load a dataset
#' @export
load_dataset <- function(datasetid, contents=list(simulation=FALSE, goldstandard=TRUE, model=TRUE, experiment=TRUE)) {
  overview <- overviewer("datasets") %>% filter(id==datasetid) %>% as.list()
  dataset <- loader(datasetid, "datasets")[[1]]
  if(contents$simulation) dataset$simulations <- loader(dataset$info$simulationid, "simulations")[[1]]
  if(contents$goldstandard) dataset$gs <- loader(dataset$info$goldstandardid, "goldstandards")[[1]]
  if(contents$model) dataset$model <- loader(dataset$info$modelid, "models")[[1]]
  if(contents$experiment) dataset <- c(dataset, loader(dataset$info$experimentid, "experiments")[[1]])
  dataset
}
