#' List all datasets
#' @export
list_datasets = function() overviewer("datasets")

#' Save a dataset of a particular type
#' @export
saver <- function(x, type, overview_only=FALSE) {
  newoverview <- map(x, ~.$info) %>% bind_rows()
  
  if(!overview_only) {
    for(xi in x) {
      path <- paste0(.datasets_location, "/", type, "/", xi$info$id, ".rds")
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(xi, path)
    }
  }
  
  overview_location <- paste0(.datasets_location, "/", type, ".rds")
  overview <- if(file.exists(overview_location)) {readRDS(paste0(.datasets_location, "/", type, ".rds"))} else {tibble()  }
  
  overview %>% 
    filter(!(id %in% newoverview$id)) %>% 
    bind_rows(newoverview) %>% 
    saveRDS(paste0(.datasets_location, "/", type, ".rds"))
}

#' Generate an overview data.frame for a particular type of data
#' @export
overviewer <- function(type) {
  overview_location <- paste0(.datasets_location, "/", type, ".rds")
  if(file.exists(overview_location)) {readRDS(paste0(.datasets_location, "/", type, ".rds"))} else {tibble(id=character())  }
}

#' Load multiple datasets of one type
#' @export
loader <- function(x, type) {
  overview <- overviewer(type)
  
  map(x, function(xi) {
    path <- paste0(.datasets_location, "/", type, "/", xi, ".rds")
    readRDS(path)
  })
}




contents_dataset <- function(simulation=FALSE, goldstandard=TRUE, model=TRUE, experiment=TRUE, special_cells=TRUE) {
  as.list(environment())
}

#' Load a dataset
#' @export
load_dataset <- function(datasetid, contents=contents_dataset()) {
  overview <- overviewer("datasets") %>% filter(id==datasetid) %>% as.list()
  dataset <- loader(datasetid, "datasets")[[1]]
  if(contents$simulation) dataset$simulations <- loader(dataset$info$simulation_id, "simulations")[[1]]
  if(contents$goldstandard) dataset$gs <- loader(dataset$info$goldstandard_id, "goldstandards")[[1]]
  if(contents$model) dataset$model <- loader(dataset$info$model_id, "models")[[1]]
  if(contents$experiment) dataset <- c(dataset, loader(dataset$info$experiment_id, "experiments")[[1]])
  
  # extract special cells for this dataset
  # we do this at dataset level because eg. the starting cell depends on which cells were selected by the platform
  if(contents$special_cells && contents$goldstandard) {
    dataset$special_cells <- list()
    
    # start cell
    # get the cell with the lowest milestone number, and within that milestone the lowest progression
    # this will not be correct in all cases, eg. multiple starting locations
    dataset$special_cells$start_cell_id <- dataset$gs$cellinfo %>% filter(cell_id %in% rownames(dataset$counts)) %>% arrange(state_id, time) %>% pull(cell_id) %>% first
    
    
    # end cell
    # get the end milestones (milestones with out-degree = 0), and select the cells with the highest percentage of those milestones
    # this will not work if there were no cells selected from the particular milestone
    # we should in the future look at the cell with the minimal geodesic distance to the end milestone
    end_milestone_ids <- dataset$gs$milestone_net$to %>% keep(~!(. %in% dataset$gs$milestone_network$from))
    dataset$special_cells$end_cell_ids <- dataset$gs$milestone_percentages %>% filter(cell_id %in% rownames(dataset$counts)) %>% filter(milestone_id %in% end_milestone_ids) %>% group_by(milestone_id) %>% summarise(cell_id = cell_id[which.max(percentage)]) %>% pull(cell_id)
  }
  
  dataset
}
