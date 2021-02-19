#' dyngen: A multi-modal simulator for spearheading single-cell omics analyses
#' 
#' A toolkit for generating synthetic single cell data. 
#' 
#' @section Step 1, initialise dyngen model:
#'   * [initialise_model()]: Define and store settings for all following steps. See each of the sections below for more information.
#'   * Use a predefined backbone: 
#'      - [list_backbones()]
#'      - [backbone_bifurcating()]
#'      - [backbone_bifurcating_converging()]
#'      - [backbone_bifurcating_cycle()]
#'      - [backbone_bifurcating_loop()]
#'      - [backbone_branching()]
#'      - [backbone_binary_tree()]
#'      - [backbone_consecutive_bifurcating()]
#'      - [backbone_trifurcating()]
#'      - [backbone_converging()]
#'      - [backbone_cycle()]
#'      - [backbone_cycle_simple()]
#'      - [backbone_linear()]
#'      - [backbone_linear_simple()]
#'      - [backbone_disconnected()]
#'   * Create a custom backbone: 
#'      - [backbone()]
#'      - [bblego()]
#'      - [bblego_linear()]
#'      - [bblego_branching()]
#'      - [bblego_start()]
#'      - [bblego_end()]
#'   * Visualise the backbone:
#'      - [plot_backbone_modulenet()]
#'      - [plot_backbone_statenet()]
#' 
#' @section Step 2, generate TF network:
#'   * [generate_tf_network()]: Generate a transcription factor network from the backbone
#'   * [tf_network_default()]: Parameters for configuring this step
#' 
#' @section Step 3, add more genes to the gene network:
#'   * [generate_feature_network()]: Generate a target network
#'   * [feature_network_default()]: Parameters for configuring this step
#'   * [plot_feature_network()]: Visualise the gene network
#' 
#' @section Step 4, generate gene kinetics:
#'   * [generate_kinetics()]: Generate the gene kinetics
#'   * [kinetics_default()]: Parameters for configuring this step
#'   * [kinetics_noise_none()], [kinetics_noise_simple()]: Different kinetics distributions to sample from
#' 
#' @section Step 5, simulate the gold standard:
#'   * [generate_gold_standard()]: Simulate the gold standard backbone, used for mapping to cell states afterwards
#'   * [gold_standard_default()]: Parameters for configuring this step
#'   * [plot_gold_mappings()]: Visualise the mapping of the simulations to the gold standard
#'   * [plot_gold_simulations()]: Visualise the gold standard simulations using the dimred
#'   * [plot_gold_expression()]: Visualise the expression of the gold standard over simulation time
#' 
#' @section Step 6, simulate the cells:
#'   * [generate_cells()]: Simulate the cells based on its GRN
#'   * [simulation_default()]: Parameters for configuring this step
#'   * [simulation_type_wild_type()], [simulation_type_knockdown()]: Used for configuring the type of simulation
#'   * [plot_simulations()]: Visualise the simulations using the dimred
#'   * [plot_simulation_expression()]: Visualise the expression of the simulations over simulation time
#' 
#' @section Step 7, simulate cell and transcripting sampling:
#'   * [generate_experiment()]: Sample cells and transcripts from experiment
#'   * [list_experiment_samplers()], [experiment_snapshot()], [experiment_synchronised()]: Parameters for configuring this step
#'   * [simtime_from_backbone()]: Determine the simulation time from the backbone
#'   * [plot_experiment_dimred()]: Plot a dimensionality reduction of the final dataset
#' 
#' @section Step 8, convert to dataset:
#'   * [as_dyno()], [wrap_dataset()]: Convert a dyngen model to a dyno dataset
#'   * [as_anndata()]: Convert da yngen model to an anndata dataset
#'   * [as_sce()]: Convert a dyngen model to a SingleCellExperiment dataset
#'   * [as_seurat()]: Convert a dyngen model to a Seurat dataset
#'   
#' @section One-shot function:
#'   * [generate_dataset()]: Run through steps 2 to 8 with a single function
#' 
#' @section Data objects:
#'   * [example_model]: A (very) small toy dyngen model, used for documentation and testing purposes
#'   * [realcounts]: A set of real single-cell expression datasets, to be used as reference datasets
#'   * [realnets]: A set of real gene regulatory networks, to be sampled in step 3
#' 
#' @section Varia functions:
#'   * [dyngen]: This help page
#'   * [get_timings()]: Extract execution timings for each of the dyngen steps
#'   * [combine_models()]: Combine multiple dyngen models
#'   * [rnorm_bounded()]: A bounded version of [rnorm()]
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
#' dataset <- as_dyno(model)
#' 
#' # dynplot::plot_dimred(dataset)
#' }
#'
#' @importFrom dplyr bind_rows filter group_by mutate mutate_all n pull sample_n select transmute ungroup do left_join last .data
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
NULL
