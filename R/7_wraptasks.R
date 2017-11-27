# Wrap tasks
#' @param params Parameters of the dataset
#' @param model Model
#' @param simulation Simulation
#' @param gs Gold standard
#' @param experiment Experiment
#' @param normalization Normalization
#' @export
wrap_task <- function(params, model, simulation, gs, experiment, normalization) {
  counts <- normalization$counts
  expression <- normalization$expression[rownames(counts), colnames(counts)]
  
  cell_ids <- rownames(counts)
  
  # create sample info
  cell_info <- experiment$cellinfo %>%
    slice(match(cell_ids, cell_id))
  
  # get milestone network
  milestone_network <- gs$milestone_network %>% 
    filter(!burn) %>% 
    mutate(length=1, directed=TRUE) %>% 
    mutate_at(c("from", "to"), as.character) %>% 
    select(from, to, length, directed)
  
  # get milestone ids
  milestone_ids <- milestone_network %>% select(from, to) %>% 
    unlist() %>% 
    unique %>%
    as.character
  
  # progressions
  progressions <- cell_info %>%
    select(step_id, cell_id) %>%
    left_join(gs$progressions, by = "step_id") %>%
    mutate_at(c("from", "to"), as.character) %>% 
    select(cell_id, from, to, percentage)
  
  # milestone percentages
  milestone_percentage <- dynutils::convert_progressions_to_milestone_percentages(cell_ids, milestone_ids, milestone_network, progressions)
  
  # add prior information
  prior_information <- dynutils::generate_prior_information(milestone_ids, milestone_network, progressions, milestone_percentage)
  
  # feature info
  feature_info <- experiment$geneinfo %>% slice(match(colnames(counts), gene_id)) %>% rename(feature_id = gene_id)
  
  # create task
  task <- dynutils::wrap_ti_task_data(
    ti_type = params$updates$modulenetname,
    id = params$updates$dataset_id,
    cell_ids = cell_ids,
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    counts = counts,
    expression = expression,
    cell_info = cell_info,
    feature_info = feature_info,
    info = params$updates
  )
  
  task$geodesic_dist <- dynutils::compute_emlike_dist(task)
  
  task
}