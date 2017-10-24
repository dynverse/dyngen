#' Load dyngen datasets meta information
#'
#' \strong{Important!} A variable \code{.datasets_location} needs to be defined in the global environment for this function
#' to work correctly.
#'
#' @return The meta information of the datasets
#' @export
#'
#' @examples
#' \dontrun{
#' .dataset_location <- "path_to_dyngen_results"
#' load_datasets_info()
#' }
load_datasets_info <- function() {
  readRDS(paste0(.datasets_location, "/datasets.rds"))
}


#' Load dyngen datasets
#'
#' \strong{Important!} A variable \code{.datasets_location} needs to be defined in the global environment for this function
#' to work correctly.
#'
#' @param datasets_info The meta information of the datasets to read
#' @param mc_cores The number of cores to use while loading in the datasets (default: 1)
#'
#' @return A tibble of datasets
#' @export
#'
#' @importFrom pbapply pblapply
#' @importFrom dynutils wrap_ti_task_data get_cell_grouping compute_emlike_dist list_as_tibble
#'
#' @examples
#' \dontrun{
#' .dataset_location <- "path_to_dyngen_results"
#' datasets <- load_datasets(mc_cores = 1, datasets_info = load_datasets_info())
#' }
load_datasets <- function(datasets_info = load_datasets_info(), mc_cores = 1) {
  # load the datasets one by one
  task_wrapped <- pbapply::pblapply(seq_len(nrow(datasets_info)), cl = mc_cores, function(dataset_num) {
    # load datasets
    dataset_id <- datasets_info$id[[dataset_num]]
    dataset <- load_dataset(dataset_id)
    
    # copy objects into environment 
    # list2env(dataset, environment())
    counts <- dataset$counts
    cellinfo <- dataset$cellinfo
    gs <- dataset$gs
    model <- dataset$model
    platform <- dataset$platform
    takesetting <- dataset$takesetting
    
    # TODO: shouldn't need this
    colnames(counts) <- paste0("Gene", seq_len(ncol(counts)))
    
    # Get the rownames
    cell_ids <- rownames(counts)
    
    # create sample info
    sample_info <- cellinfo %>%
      slice(match(cell_ids, step_id)) %>%
      select(cell_id, step, simulation_time = simulationtime)
    
    # create milestones
    milestone_ids <- gs$milestone_percentages$milestone_id %>%
      unique %>%
      as.character
    
    # create network
    milestone_network <- gs$milestone_network %>%
      mutate_at(c("from", "to"), as.character)
    
    # create progressions
    progression <- sample_info %>%
      select(cell_id) %>%
      left_join(gs$progression, by = "cell_id") %>%
      mutate_at(c("from", "to"), as.character)
    
    # create task
    task <- dynutils::wrap_ti_task_data(
      ti_type = model$modulenetname,
      id = dataset_id,
      cell_ids = cell_ids,
      milestone_ids = milestone_ids,
      milestone_network = milestone_network,
      progressions = progression,
      counts = counts,
      sample_info = sample_info,
      task_ix = dataset_num,
      modulenet_id = model$modulenetname,
      platform_id = platform$platform_id,
      takesetting_type = takesetting$type,
      model_replicate = model$modelsetting$replicate,
      special_cells = special_cells
    )
    
    # add cell grouping
    task$cell_grouping <- dynutils::get_cell_grouping(task$milestone_percentages)
    
    # add geodesic dist
    task$geodesic_dist <- dynutils::compute_emlike_dist(task)
    
    # return task
    task
  })
  
  # return tasks
  dynutils::list_as_tibble(task_wrapped) %>%
    left_join(datasets_info, by = c("id" = "id"))
}
