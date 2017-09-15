#' Convert milestone percentages to progressions
#' @param cell_ids Vector of all cell ids
#' @param milestone_ids Vector of milestone ids
#' @param milestone_network Milestone network
#' @param milestone_percentages Milestone percentages
#' @value Progressions
#' @export
convert_milestone_percentages_to_progressions <- function(cell_ids, milestone_ids, milestone_network, milestone_percentages) {
  bind_rows(lapply(cell_ids, function(cid) {
    relevant_pct <- milestone_percentages %>% filter(cell_id == cid)

    if (nrow(relevant_pct) >= 2) {
      relevant_progr <- milestone_network %>%
        filter(from %in% relevant_pct$milestone_id & to %in% relevant_pct$milestone_id) %>%
        left_join(relevant_pct, by = c("to" = "milestone_id")) %>%
        select(cell_id, from, to, percentage)
      if (nrow(relevant_progr) == 0) {
        stop("According to milestone_percentages, cell ", sQuote(cid), " is between milestones ",
             paste(sQuote(relevant_pct$milestone_id), collapse = " and "), ", but this edge does not exist in milestone_network!")
      }
    } else if (nrow(relevant_pct) == 1) {
      relevant_net <- milestone_network %>% filter(to %in% relevant_pct$milestone_id)
      if (nrow(relevant_net) == 0) {
        relevant_net <- milestone_network %>% filter(from %in% relevant_pct$milestone_id)
      }
      relevant_progr <- relevant_net %>%
        mutate(cell_id = cid, percentage = 1) %>%
        select(cell_id, from, to, percentage)
    } else {
      relevant_progr <- NULL
    }
    relevant_progr
  }))

}


#' Convert progressions to milestone percentages
#' @param cell_ids Vector of all cell ids
#' @param milestone_ids Vector of milestone ids
#' @param milestone_network Milestone network
#' @param progressions Progressions dataframe
#' @value Milestone percentages
#' @export
convert_progressions_to_milestone_percentages <- function(cell_ids, milestone_ids, milestone_network, progressions) {
  check_froms <- progressions %>% group_by(cell_id) %>% summarise(n = length(unique(from)))
  if (any(check_froms$n > 1)) {
    stop("In ", sQuote("progressions"), ", cells should only have 1 unique from milestone.")
  }

  froms <- progressions %>% group_by(cell_id) %>% summarise(milestone_id = from[[1]], percentage = 1 - sum(percentage))
  tos <- progressions %>% select(cell_id, milestone_id = to, percentage)
  bind_rows(froms, tos)
}

add_phantom_edges <- function(milestone_ids, milestone_network) {
  bind_rows(
    milestone_network,
    bind_rows(lapply(milestone_ids, function(x) {
      strx <- milestone_network %>%
        filter(from == x)

      if (nrow(strx) > 1) {
        strx <- strx %>%
          mutate(
            angle = seq(0, 120/360*pi*2, length.out = n()),
            x = length * cos(angle),
            y = length * sin(angle)
          )
        poss <- strx %>% select(x, y) %>% as.matrix
        rownames(poss) <- strx$to
        poss %>%
          dist %>%
          as.matrix %>%
          reshape2::melt(varnames = c("from", "to"), value.name = "length") %>%
          mutate(from = as.character(from), to = as.character(to)) %>%
          filter(from != to)
      } else {
        NULL
      }
    })),
    bind_rows(lapply(milestone_ids, function(x) {
      strx <- milestone_network %>%
        filter(to == x)

      if (nrow(strx) > 1) {
        strx <- strx %>%
          mutate(
            angle = seq(0, 120/360*pi*2, length.out = n()),
            x = length * cos(angle),
            y = length * sin(angle)
          )
        poss <- strx %>% select(x, y) %>% as.matrix
        rownames(poss) <- strx$from
        poss %>%
          dist %>%
          as.matrix %>%
          reshape2::melt(varnames = c("from", "to"), value.name = "length") %>%
          mutate(from = as.character(from), to = as.character(to)) %>%
          filter(from != to)
      } else {
        NULL
      }
    })))
}

