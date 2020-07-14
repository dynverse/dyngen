#' dyngen: Synthetic single cell data
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
#' @importFrom furrr future_map
#' @importFrom rlang %|%
#'
#' @docType package
#' @name dyngen
NULL