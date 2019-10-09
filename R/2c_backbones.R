#' @export
#' @rdname backbone_models
backbone_bifurcating <- function() {
  bblego(
    bblego_start("A", type = "simple"),
    bblego_linear("A", "B"),
    bblego_branching("B", c("C", "D")),
    bblego_end("C"),
    bblego_end("D")
  )
}


#' @export
#' @rdname backbone_models
backbone_bifurcating_converging <- function() {
  module_info <- tribble(
    ~module_id, ~ba, ~burn,
    "A1", 1, TRUE,
    "B1", 0, FALSE,
    "B2", 0, FALSE,
    "C1", 0, FALSE,
    "C2", 0, FALSE,
    "C3", 0, FALSE,
    "C4", 0, FALSE,
    "D1", 0, FALSE,
    "D2", 0, FALSE,
    "D3", 0, FALSE,
    "D4", 0, FALSE,
    "E1", 0, FALSE,
    "F1", 0, FALSE,
    "F2", 0, FALSE
  )
  
  module_network <- bind_rows(
    modnet_chain(c("A1", "B1", "B2")),
    modnet_edge("B2", c("C1", "D1")),
    modnet_pairwise(c("C1", "D1"), effect = -1, strength = 100),
    modnet_edge(c("C1", "D1"), c("C2", "D2"), strength = 2),
    modnet_edge(c("C2", "D2"), c("C1", "D1"), strength = 2),
    modnet_chain(c("C2", "C3", "C4", "E1")),
    modnet_chain(c("D2", "D3", "D4", "E1")),
    modnet_pairwise(c("C1", "D1"), c("C2", "D2"), effect = -1, strength = 100000),
    modnet_edge("E1", c("C1", "C2", "D1", "D2"), effect = -1, strength = 100),
    modnet_edge("E1", c("C4", "D4")),
    modnet_chain(c("E1", "F1", "F2"))
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "sBurn", "sA", "+A1", TRUE, TRUE, 2,
    "sA", "sB", "+B1,+B2", FALSE, FALSE, 2,
    "sB", "sC", "+C1,+C2,+C3,+C4", FALSE, FALSE, 2,
    "sB", "sD", "+D1,+D2,+D3,+D4", FALSE, FALSE, 2,
    "sC", "sE", "+E1,+D1,+D2,+D3,+D4", FALSE, FALSE, 4,
    "sD", "sE", "+E1,+C1,+C3,+C3,+C4", FALSE, FALSE, 4,
    "sE", "sF", "+F1,+F2", FALSE, FALSE, 2
  )
  
  backbone(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname backbone_models
backbone_bifurcating_cycle <- function() {
  module_info <- tribble(
    ~module_id, ~ba, ~burn,
    "A1", 1, TRUE,
    "B1", 0, FALSE,
    "B2", 0, FALSE,
    "A2", 1, TRUE,
    "A3", 1, TRUE,
    "C1", 0, TRUE,
    "C2", 0, TRUE,
    "C3", 0, TRUE,
    "D1", 0, TRUE,
    "D2", 0, TRUE,
    "D3", 0, TRUE,
    "E1", 0, TRUE,
    "E2", 0, TRUE,
    "E3", 0, TRUE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "A1", "B1", 1, 1, 2,
    "B1", "B2", 1, 1, 2,
    "B1", "A2", -1, 10, 2,
    "B2", "A3", -1, 10, 2,
    "B2", "C1", 1, 1, 2, 
    "B2", "D1", 1, 1, 2,
    "C1", "C1", 1, 4, 2,
    "D1", "D1", 1, 4, 2,
    "C1", "D1", -1, 100, 2,
    "D1", "C1", -1, 100, 2,
    "C1", "C2", 1, 1, 2,
    "C2", "C3", 1, 1, 2,
    "D1", "D2", 1, 1, 2,
    "D2", "D3", 1, 1, 2,
    "C3", "E1", 1, 1, 2,
    "D3", "E1", 1, 1, 2,
    "E1", "E1", 1, 5, 2,
    "E1", "A1", -1, 10, 2,
    "E1", "C1", -1, 20, 2,
    "E1", "D1", -1, 20, 2,
    "E1", "E2", 1, 1, 2,
    "E2", "E3", 1, 1, 2,
    "E3", "E1", -1, 10, 2,
    "E3", "C1", -1, 20, 2,
    "E3", "D1", -1, 20, 2,
    "E3", "A1", 1, 2, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "sBurn", "sB", "+A1,+A2,+A3,+B1,+B2", TRUE, TRUE, 3,
    "sB", "sC", "+C1,+C2,+C3", FALSE, FALSE, 2,
    "sB", "sD", "+D1,+D2,+D3", FALSE, FALSE, 2,
    "sC", "sE", "+E1", FALSE, FALSE, 4,
    "sD", "sE", "+E1", FALSE, FALSE, 4,
    "sE", "sA", "+E2,+E3,-B1,-B2,-C1,-C2,-C3,-D1,-D2,-D3", FALSE, FALSE, 4,
    "sA", "sB", "+B1,+B2,-E1,-E2|-E3", FALSE, FALSE, 2
  )
  
  backbone(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname backbone_models
backbone_bifurcating_loop <- function() {
  module_info <- tribble(
    ~module_id, ~ba, ~burn,
    "A1", 1, TRUE,
    "A2", 0, TRUE,
    "A3", 1, TRUE,
    "B1", 0, FALSE,
    "B2", 1, TRUE,
    "C1", 0, FALSE,
    "C2", 0, FALSE,
    "C3", 0, FALSE,
    "D1", 0, FALSE,
    "D2", 0, FALSE,
    "D3", 1, TRUE,
    "D4", 0, FALSE,
    "D5", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "A1", "A2", 1, 10, 2,
    "A2", "A3", -1,  10, 2,
    "A2", "B1", 1, 1, 2,
    "B1", "B2", -1, 10, 2,
    "B1", "C1", 1, 1, 2,
    "B1", "D1", 1, 1, 2,
    "C1", "C1", 1, 10, 2,
    "C1", "D1", -1, 100, 2,
    "C1", "C2", 1, 1, 2,
    "C2", "C3", 1, 1, 2,
    "C2", "A2", -1, 10, 2,
    "C2", "B1", -1, 10, 2,
    "C3", "A1", -1, 10, 2,
    "C3", "C1", -1, 10, 2,
    "C3", "D1", -1, 10, 2,
    "D1", "D1", 1, 10, 2,
    "D1", "C1", -1, 100, 2,
    "D1", "D2", 1, 1, 2,
    "D1", "D3", -1, 10, 2,
    "D2", "D4", 1, 1, 2,
    "D4", "D5", 1, 1, 2,
    "D3", "D5", -1, 10, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "sBurn", "sA", "+A1,+A2,+A3,+B2,+D3", TRUE, TRUE, 2,
    "sA", "sB", "+B1", FALSE, FALSE, 2,
    "sB", "sC", "+C1,+C2|-A2,-B1,+C3|-C1,-D1,-D2", FALSE, FALSE, 3,
    "sB", "sD", "+D1,+D2,+D4,+D5", FALSE, FALSE, 4,
    "sC", "sA", "+A1,+A2", FALSE, FALSE, 2
  )
  
  backbone(module_info, module_network, expression_patterns)
}

#' @param num_modifications The number of branch points in the generated backbone.
#' @param min_degree The minimum degree of each node in the backbone.
#' @param max_degree The maximum degree of each node in the backbone.
#' 
#' @export
#' @rdname backbone_models
#' @importFrom stats rbinom
backbone_branching <- function(
  num_modifications = rbinom(1, size = 6, 0.25) + 1,
  min_degree = 3,
  max_degree = sample(min_degree:5, 1)
) {
  assert_that(
    num_modifications >= 0,
    min_degree >= 3,
    max_degree >= 3,
    min_degree <= max_degree
  )
  
  statenet <- tribble(
    ~from, ~to,
    "A", "B",
    "B", "C"
  )
  num_nodes <- 3
  
  for (i in seq_len(num_modifications)) {
    j <- if (nrow(statenet) == 2) 2 else sample(2:nrow(statenet), 1)
    from_ <- statenet$from[[j]]
    to_ <- statenet$to[[j]]
    num_new_nodes <- rbinom(1, size = max_degree - min_degree, 0.25) + min_degree - 1
    new_nodes <- LETTERS[num_nodes + seq_len(num_new_nodes)]
    num_nodes <- num_nodes + num_new_nodes
    statenet <-
      bind_rows(
        statenet %>% slice(seq_len(j-1)),
        tibble(from = from_, to = new_nodes[[1]]),
        statenet %>% slice(j) %>% mutate(from = new_nodes[[1]]),
        tibble(from = new_nodes[[1]], to = new_nodes[-1]),
        if (j != nrow(statenet)) statenet %>% slice(seq(j+1, nrow(statenet))) else NULL
      )
  }
  
  state_names <- statenet %>% as.matrix %>% t %>% as.vector %>% unique()
  
  state_name_map <- set_names(sort(state_names), state_names)
  statenet <- statenet %>% 
    mutate(
      from = ifelse(from %in% names(state_name_map), state_name_map[from], from),
      to = ifelse(to %in% names(state_name_map), state_name_map[to], to)
    )
  
  sngr <- 
    statenet %>%
    group_by(from) %>% 
    summarise(tos = list(to)) %>% 
    mutate(from = factor(from, state_name_map)) %>% 
    arrange(from)
    
  
  bblego_list <- c(
    list(bblego_start(statenet %>% filter(!from %in% to) %>% pull(from), type = "simple")),
    map(
      seq_len(nrow(sngr)), 
      function(i) {
        from_node <- sngr$from[[i]] %>% as.character()
        to_node <- sngr$tos[[i]]
        
        if (length(to_node) == 1) {
          bblego_linear(from = from_node, to = to_node)
        } else {
          bblego_branching(from = from_node, to = to_node)
        }
      }
    ),
    map(
      statenet %>% filter(!to %in% from) %>% pull(to) %>% sort(), 
      bblego_end
    )
  )
  
  bblego(.list = bblego_list)
}

#' @export
#' @rdname backbone_models
#' @importFrom stats rbinom
backbone_binary_tree <- function(
  num_modifications = rbinom(1, size = 6, 0.25) + 1
) {
  backbone_branching(num_modifications = num_modifications, min_degree = 3, max_degree = 3)
}


#' @export
#' @rdname backbone_models
backbone_consecutive_bifurcating <- function() {
  backbone_branching(num_modifications = 2, max_degree = 3)
}

#' @export
#' @rdname backbone_models
backbone_trifurcating <- function() {
  backbone_branching(num_modifications = 1, min_degree = 4, max_degree = 4)
}

#' @export
#' @rdname backbone_models
backbone_converging <- function() {
  module_info <- tribble(
    ~module_id, ~ba, ~burn,
    "A1", 1, TRUE,
    "A2", 1, TRUE,
    "A3", 1, TRUE,
    "B1", 0, TRUE,
    "B2", 0, TRUE,
    "B3", 0, FALSE,
    "C1", 0, TRUE,
    "C2", 0, TRUE,
    "C3", 0, FALSE,
    "D1", 0, FALSE,
    "D2", 1, TRUE,
    "E1", 1, FALSE,
    "E2", 0, FALSE,
    "E3", 0, FALSE,
    "E4", 0, FALSE,
    "E5", 0, FALSE
  )
  
  module_network <- bind_rows(
    modnet_chain(c("A1", "A2", "A3")),
    modnet_edge("A3", c("B1", "C1")),
    modnet_pairwise(c("B1", "C1"), effect = -1, strength = 100),
    modnet_edge(c("B1", "C1"), c("B2", "C2"), strength = 2),
    modnet_edge(c("B2", "C2"), c("B1", "C1"), strength = 2),
    modnet_chain(c("B2", "B3", "D1")),
    modnet_chain(c("C2", "C3", "D1")),
    modnet_edge(c("D1", "D1"), c("B3", "C3")),
    modnet_pairwise(c("B1", "C1"), c("B2", "C2"), effect = -1, strength = 100000),
    modnet_edge("D1", c("B1", "C1", "B2", "C2"), effect = -1, strength = 100),
    modnet_edge(c("B2", "C2"), "E1", effect = -1, strength = 1),
    modnet_chain(c("D1", "D2", "E1"), effect = -1, strength = 1),
    modnet_chain(c("E1", "E2", "E3", "E4", "E5"))
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "sBurn1", "sBurn2", "+A1,+A2,+A3,+D2", TRUE, TRUE, 5,
    "sBurn2", "preB", "+B1,+B2", FALSE, TRUE, 3,
    "sBurn2", "preC", "+C1,+C2", FALSE, TRUE, 3,
    "preB", "sB", "+B3", FALSE, FALSE, 3,
    "preC", "sC", "+C3", FALSE, FALSE, 3,
    "sB", "sD", "+D1,+C1,+C2,+C3", FALSE, FALSE, 3,
    "sC", "sD", "+D1,+B1,+B2,+B3", FALSE, FALSE, 3,
    "sD", "sE", "+E1,+E2,+E3", FALSE, FALSE, 3
  )
  
  backbone(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname backbone_models
backbone_cycle <- function() {
  module_info <- tribble(
    ~module_id, ~ba, ~burn,
    "M1", 1, TRUE,
    "M2", 1, TRUE,
    "M3", 1, TRUE,
    "M4", 0, FALSE,
    "M5", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "M1", "M2", -1, 4, 2,
    "M2", "M3", -1, 4, 2,
    "M3", "M4", 1, 1, 2,
    "M4", "M5", 1, 1, 2,
    "M5", "M1", -1, 4, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "sBurn", "s1", "+M1,+M2,+M3", TRUE, TRUE, 3,
    "s1", "s2", "+M4,+M5,-M3", FALSE, FALSE, 4,
    "s2", "s3", "+M3,-M1,-M2", FALSE, FALSE, 4,
    "s3", "s1", "+M1,+M2,-M4,-M5", FALSE, FALSE, 4
  )
  
  backbone(module_info, module_network, expression_patterns)
}

#' @export
#' @rdname backbone_models
backbone_linear <- function() {
  bblego(
    bblego_start("A", type = "simple"),
    bblego_linear("A", "B"),
    bblego_linear("B", "C"),
    bblego_end("C")
  )
}



#' @export
#' @rdname backbone_models
#' @param left_backbone A backbone (other than a disconnected backbone), see [list_backbones()].
#' @param right_backbone A backbone (other than a disconnected backbone), see [list_backbones()].
#' @param num_common_modules The number of modules which are regulated by either backbone.
#' 
#' @importFrom methods is
backbone_disconnected <- function(
  left_backbone = list_backbones() %>% keep(., names(.) != "disconnected") %>% sample(1) %>% first(),
  right_backbone = list_backbones() %>% keep(., names(.) != "disconnected") %>% sample(1) %>% first(),
  num_common_modules = 10
) {
  if (is.character(left_backbone)) {
    left_backbone <- list_backbones()[[left_backbone]]
  }
  if (is.character(right_backbone)) {
    right_backbone <- list_backbones()[[right_backbone]]
  }
  if (is.function(left_backbone)) {
    left_backbone <- left_backbone()
  }
  if (is.function(right_backbone)) {
    right_backbone <- right_backbone()
  }
  assert_that(
    is(left_backbone, "dyngen::backbone"),
    is(right_backbone, "dyngen::backbone")
  )
  
  left <- left_backbone
  right <- right_backbone
  
  leftpas <- function(x) paste0("left_", x)
  rightpas <- function(x) paste0("right_", x)
  
  lmi <- left$module_info %>% mutate_at(vars(module_id), leftpas)
  rmi <- right$module_info %>% mutate_at(vars(module_id), rightpas)
  
  lmn <- left$module_network %>% mutate_at(vars(from, to), leftpas)
  rmn <- right$module_network %>% mutate_at(vars(from, to), rightpas)
  
  lep <- left$expression_patterns %>% mutate_at(vars(from, to), leftpas) %>% mutate(module_progression = gsub("([+-])", "\\1left_", module_progression))
  rep <- right$expression_patterns %>% mutate_at(vars(from, to), rightpas) %>% mutate(module_progression = gsub("([+-])", "\\1right_", module_progression))
  
  # find starting modules and states; 
  # assume there is only one starting module / state per backbone
  lmis <- lmi %>% filter(ba > 0) %>% pull(module_id)
  rmis <- rmi %>% filter(ba > 0) %>% pull(module_id)
  leps <- lep %>% filter(start) %>% pull(from)
  reps <- rep %>% filter(start) %>% pull(from)
  
  assert_that(
    length(lmis) > 0,
    length(rmis) > 0,
    length(leps) == 1,
    length(reps) == 1
  )
  
  common_modules <- paste0("B", seq_len(num_common_modules))
  
  module_info <- bind_rows(
    tibble(
      module_id = paste0("A", 1:7),
      ba = ifelse(module_id == "A1", 1, 0),
      burn = TRUE,
      ind = 1
    ),
    lmi %>% mutate(ba = 0),
    rmi %>% mutate(ba = 0),
    tibble(module_id = common_modules, ba = 0, burn = FALSE)
  ) %>% select(-color)
  
  module_network <- bind_rows(
    modnet_chain(c("A1", "A2", "A3")),
    modnet_edge("A3", c("A4", "A5")),
    modnet_pairwise(c("A4", "A5"), effect = -1, strength = 100),
    modnet_edge(c("A4", "A5", "A6", "A7"), c("A6", "A7", "A4", "A5"), strength = 2),
    modnet_pairwise(c("A4", "A5"), c("A6", "A7"), effect = -1, strength = 100000),
    modnet_edge(c("A6", "A7"), "A1", effect = -1, strength = 10),
    lmi %>% filter(ba > 0) %>% transmute(from = "A6", to = module_id, effect = 1, strength = ba, cooperativity = 2),
    rmi %>% filter(ba > 0) %>% transmute(from = "A7", to = module_id, effect = 1, strength = ba, cooperativity = 2),
    lmn,
    rmn,
    tibble(
      from = c(
        sample(lmn$from, length(common_modules), replace = TRUE),
        sample(rmn$from, length(common_modules), replace = TRUE)
      ),
      to = c(common_modules, common_modules),
      effect = 1,
      strength = 1,
      cooperativity = 2
    )
  )
  
  expression_patterns <-
    bind_rows(
      tribble(
        ~from, ~to, ~module_progression, ~start, ~burn, ~time,
        "sBurn", "sA", paste0("+", c("A1", "A2", "A3", common_modules), collapse = ","), TRUE, TRUE, 2,
        "sA", leps, "+A4,+A6", FALSE, TRUE, 2,
        "sA", reps, "+A5,+A7", FALSE, TRUE, 2
      ),
      lep %>% mutate(start = FALSE),
      rep %>% mutate(start = FALSE)
    )
  
  backbone(module_info, module_network, expression_patterns)
}
