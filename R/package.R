#' dyngen: A multi-modal simulator for spearheading single-cell omics analyses
#' 
#' A toolkit for generating synthetic single cell data. 
#'
#' @importFrom dplyr bind_rows filter group_by mutate mutate_all n pull sample_n select transmute ungroup do left_join last
#' @importFrom dplyr row_number bind_cols full_join summarise inner_join slice rename case_when arrange first mutate_at vars nth
#' @importFrom tidyr gather unnest everything one_of crossing
#' @importFrom purrr %>% map map_df map_chr keep pmap map2 set_names map_int map_dbl
#' @importFrom tibble as_tibble tibble enframe deframe lst tribble
#' @importFrom methods as
#' @importFrom utils head tail
#' @importFrom assertthat assert_that %has_name%
#' @importFrom dynutils %all_in% extract_row_to_list is_sparse %has_names% add_class extend_with
#' @importFrom Matrix t Matrix sparseMatrix summary
#' @importFrom pbapply pblapply 
#' @importFrom rlang %|% %||%
#'
#' @docType package
#' @name dyngen
#' 
#' @examples 
#' model <- initialise_model(
#'   backbone = backbone_bifurcating()
#' )
#' \dontshow{
#' # actually use a smaller example 
#' # to reduce execution time during
#' # testing of the examples
#' model <- initialise_model(
#'   backbone = model$backbone,
#'   num_cells = 5,
#'   num_targets = 0,
#'   num_hks = 0,
#'   gold_standard_params = gold_standard_default(census_interval = 1, tau = 0.1),
#'   simulation_params = simulation_default(
#'     burn_time = 10,
#'     total_time = 10,
#'     census_interval = 1,
#'     ssa_algorithm = ssa_etl(tau = 0.1),
#'     experiment_params = simulation_type_wild_type(num_simulations = 1)
#'   )
#' )
#' }
#' \donttest{
#' model <- model %>%
#'   generate_tf_network() %>%
#'   generate_feature_network() %>%
#'   generate_kinetics() %>%
#'   generate_gold_standard() %>%
#'   generate_cells() %>%
#'   generate_experiment()
#'   
#' dataset <- wrap_dataset(model)
#' 
#' # dynplot::plot_dimred(dataset)
#' }
NULL