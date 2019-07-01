#' Generate a dataset
#' 
#' This function contains the complete pipeline for generating a dataset
#' with \pkg{dyngen}. In order to have more control over how the dataset 
#' is generated, run each of the steps in this function separately.
#' 
#' @param model A dyngen initial model created with [initialise_model()].
#' 
#' @export
generate_dataset <- function(model) {
  model <- model %>% 
    generate_tf_network() %>% 
    generate_feature_network() %>% 
    generate_kinetics() %>% 
    generate_gold_standard() %>%  
    generate_cells() %>% 
    generate_experiment()
  
  traj <- 
    model %>% 
    wrap_dataset()
  
  lst(
    model,
    traj
  )
}
