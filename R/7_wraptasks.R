#' Wrap tasks
#' @param params Parameters of the dataset
#' @param model Model
#' @param simulation Simulation
#' @param gs Gold standard
#' @param experiment Experiment
#' @param normalisation Normalisation
#' 
#' @importFrom dynnormaliser generate_prior_information
#' 
#' @export
wrap_task <- function(id = "", params, model, simulation, gs, experiment, normalisation) {
  counts <- normalisation$counts
  expression <- normalisation$expression[rownames(counts), colnames(counts)]
  
  cell_ids <- rownames(counts)
  
  # create sample info
  cell_info <- experiment$cellinfo %>%
    slice(match(cell_ids, cell_id)) %>% 
    left_join(simulation$stepinfo, by="step_id")
  
  # get milestone network
  milestone_network <- gs$milestone_network %>% 
    filter(!burn) %>% 
    mutate(length=1, directed=TRUE) %>% 
    mutate_at(c("from", "to"), as.character) %>% 
    select(from, to, length, directed)
  
  # progressions
  progressions <- cell_info %>%
    select(step_id, cell_id) %>%
    left_join(gs$progressions, by = "step_id") %>%
    mutate_at(c("from", "to"), as.character) %>% 
    select(cell_id, from, to, percentage)
  
  # filter milestone_network for those edges present in the data
  milestone_network <- progressions %>% 
    select(from, to) %>% 
    distinct(from, to) %>% 
    left_join(milestone_network, c("from", "to"))
  
  # get milestone ids
  milestone_ids <- milestone_network %>% select(from, to) %>% 
    unlist() %>% 
    unique %>%
    as.character
  
  # milestone percentages
  milestone_percentage <- convert_progressions_to_milestone_percentages(cell_ids, milestone_ids, milestone_network, progressions)
  
  # feature info
  feature_info <- experiment$geneinfo %>% slice(match(colnames(counts), gene_id))
  feature_info$feature_id <- feature_info$gene_id
  
  # create task
  wrap_data(
    id = id,
    cell_ids = cell_ids,
    cell_info = cell_info,
    task_source = "synthetic",
    settings = params$settings,
    normalisation_info = normalisation$info
  ) %>% add_trajectory_to_wrapper(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    divergence_regions = NULL, # dyngen does not support divergence regions right now
    progressions = progressions
  ) %>% add_expression_to_wrapper(
    counts = counts,
    expression = expression,
    feature_info = feature_info
  ) %>% dynnormaliser::add_prior_information_to_wrapper()
}
