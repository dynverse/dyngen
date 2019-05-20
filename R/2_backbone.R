#' Backbone of the simulation model
#' 
#' A module is a group of genes which, to some extent, shows the same
#' expression behaviour. Several modules are connected together such that
#' one or more genes from one module will regulate the expression of
#' another module. By creating chains of modules, a dynamic behaviour in gene 
#' regulation can be created.
#' 
#' A backbone model contains three tibbles, namely `module_info`,
#' `module_network`, and `expression_patterns`.
#' 
#' `module_info` contains meta information on the modules themselves.
#' 
#' * module_id (character): the name of the module
#' * a0 (numeric): basal expression level of genes in this module
#' * burn (logical): whether or not outgoing edges of this module will 
#'   be active during the burn in phase
#'     
#' `module_network` describes which modules regulate which other modules.
#' 
#' * from (character): the regulating module
#' * to (character): the target module
#' * effect (integer): `1` if the regulating module upregulates 
#'   the target module, `-1` if it downregulates
#' * strength (numeric): the strength of the interaction
#' * cooperativity (numeric): cooperativity factor, larger 1 if positive cooperativity,
#'   between 0 and 1 for negative cooperativity
#' 
#' `expression_patterns` describes the expected expression pattern
#' changes when a cell is simulated by dyngen.
#' 
#' * from (character): name of a cell state
#' * to (character): name of a cell state
#' * module_progression (character): differences in module expression between the two states. 
#'   Example: `"+4,-1|+9|-4"` means the expression of module 4 will go up at the same time as module 1 goes down; 
#'   afterwards module 9 expression will go up, and afterwards module 4 expression will go down again.
#' * start (logical): Whether or not this from cell state is the start of the trajectory
#' * burn (logical): Whether these cell states are part of the burn in phase. Cells will 
#'   not get sampled from these cell states.
#' 
#' @importFrom grDevices rainbow
#' @export
backbone <- function(
  module_info,
  module_network,
  expression_patterns
) {
  assert_that(
    module_network$from %all_in% module_info$module_id,
    module_network$to %all_in% module_info$module_id
  )
  
  tmp_modules <-
    expression_patterns$module_progression %>%
    strsplit("[,|]") %>%
    unlist() %>%
    unique() %>% 
    gsub("[+-]", "", .)
  assert_that(tmp_modules %all_in% module_info$module_id)
  
  if (! module_info %has_name% "color") {
    module_info <-
      module_info %>% 
      mutate(color = grDevices::rainbow(n()))
  }
  
  lst(
    module_info, 
    module_network, 
    expression_patterns
  ) %>% 
    add_class("dyngen::backbone")
}

#' @export
#' @rdname backbone
list_backbones <- function() {
  list(
    bifurcating = backbone_bifurcating,
    bifurcating_converging = backbone_bifurcating_converging,
    bifurcating_cycle = backbone_bifurcating_cycle,
    bifurcating_loop = backbone_bifurcating_loop,
    binary_tree = backbone_binary_tree,
    branching = backbone_branching,
    consecutive_bifurcating = backbone_consecutive_bifurcating,
    converging = backbone_converging,
    cycle = backbone_cycle,
    disconnected = backbone_disconnected,
    linear = backbone_linear,
    trifurcating = backbone_trifurcating
  )
}

#' @export
#' @rdname backbone
backbone_bifurcating <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "A1", 1, TRUE,
    "B1", 0, FALSE,
    "B2", 0, FALSE,
    "C1", 0, FALSE,
    "C2", 0, FALSE,
    "D1", 0, FALSE,
    "D2", 1, TRUE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "A1", "B1", 1, 1, 2,
    "B1", "B2", 1, 1, 2,
    "B2", "C1", 1, 1, 2,
    "B2", "D1", 1, 1, 2,
    "C1", "C1", 1, 2, 2,
    "D1", "D1", 1, 2, 2,
    "C1", "D1", -1, 100, 2,
    "D1", "C1", -1, 100, 2,
    "C1", "C2", 1, 1, 2,
    "D1", "D2", -1, 3, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "sBurn", "sA", "+A1,+D2", TRUE, TRUE, 2,
    "sA", "sB", "-D2,+B1,+B2", FALSE, FALSE, 1,
    "sB", "sC", "+C1,+C2", FALSE, FALSE, 3,
    "sB", "sD", "+D1,+D2", FALSE, FALSE, 3
  )
  
  backbone(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname backbone
backbone_bifurcating_converging <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "A1", 1, TRUE,
    "B1", 0, FALSE,
    "B2", 0, FALSE,
    "C1", 0, FALSE,
    "C2", 0, FALSE,
    "C3", 0, FALSE,
    "D1", 0, FALSE,
    "D2", 0, FALSE,
    "D3", 0, FALSE,
    "E1", 0, FALSE,
    "F1", 0, FALSE,
    "F2", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "A1", "B1", 1, 1, 2,
    "B1", "B2", 1, 1, 2,
    "B2", "C1", 1, 1, 2,
    "B2", "D1", 1, 1, 2,
    "C1", "C1", 1, 10, 2,
    "C1", "C2", 1, 1, 5,
    "C1", "D1", -1, 100, 2,
    "C2", "C3", 1, 1, 5,
    "C3", "E1", 1, 1, 2,
    "D1", "D1", 1, 10, 2,
    "D1", "C1", -1, 100, 2,
    "D1", "D2", 1, 1, 5,
    "D2", "D3", 1, 1, 5,
    "D3", "E1", 1, 1, 2,
    "E1", "E1", 1, 10, 2,
    "E1", "C1", -1, 100, 2,
    "E1", "D1", -1, 100, 2,
    "E1", "F1", 1, 1, 2,
    "F1", "B2", -1, 10, 2,
    "F1", "F2", 1, 1, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "sBurn", "sA", "+A1", TRUE, TRUE, 2,
    "sA", "sB", "+B1,+B2", FALSE, FALSE, 2,
    "sB", "sC", "+C1,+C2,+C3", FALSE, FALSE, 2,
    "sB", "sD", "+D1,+D2,+D3", FALSE, FALSE, 2,
    "sC", "sE", "+E1", FALSE, FALSE, 4,
    "sD", "sE", "+E1", FALSE, FALSE, 4,
    "sE", "sF", "+F1,+F2", FALSE, FALSE, 2
  )
  
  backbone(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname backbone
backbone_bifurcating_cycle <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
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
    "C1", "C1", 1, 10, 2,
    "D1", "D1", 1, 10, 2,
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
    "sBurn", "sB", "+A1,+A2,+A3,+B1,+B2", TRUE, TRUE, 2,
    "sB", "sC", "+C1,+C2,+C3", FALSE, FALSE, 2,
    "sB", "sD", "+D1,+D2,+D3", FALSE, FALSE, 2,
    "sC", "sE", "+E1", FALSE, FALSE, 4,
    "sD", "sE", "+E1", FALSE, FALSE, 4,
    "sE", "sA", "+E2,+E3,-B1,-B2,-C1,-C2,-C3,-D1,-D2,-D3", FALSE, FALSE, 4,
    "sA", "sB", "+B1,+B2", FALSE, FALSE, 2
  )
  
  backbone(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname backbone
backbone_bifurcating_loop <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
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

#' @export
#' @rdname backbone
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

  milnet <- tribble(
    ~from, ~to, ~modules,
    "sBurn", "sA", paste0("A", seq_len(3)),
    "sA", "sB", paste0("B", seq_len(2))
  )
  num_nodes <- 2

  for (i in seq_len(num_modifications)) {
    j <- if (nrow(milnet) == 2) 2 else sample(2:nrow(milnet), 1)
    from_ <- milnet$from[[j]]
    to_ <- milnet$to[[j]]
    num_new_nodes <- rbinom(1, size = max_degree - min_degree, 0.25) + min_degree - 1
    new_nodes <- paste0("s", LETTERS[num_nodes + seq_len(num_new_nodes)])
    num_nodes <- num_nodes + num_new_nodes
    milnet <-
      bind_rows(
        milnet %>% slice(seq_len(j-1)),
        tibble(from = from_, to = new_nodes[[1]], modules = list(paste0(gsub("s", "", new_nodes[[1]]), seq_len(1)))),
        milnet %>% slice(j) %>% mutate(from = new_nodes[[1]]),
        tibble(from = new_nodes[[1]], to = new_nodes[-1], modules = map(to, ~ paste0(gsub("s", "", .), seq_len(sample(2:3, 1, prob = c(.75, .25)))))),
        if (j != nrow(milnet)) milnet %>% slice(seq(j+1, nrow(milnet))) else NULL
      )
  }
  
  # determine which modules are repressing
  milnet <- milnet %>% 
    mutate(
      double_repr =
        from != "sBurn" &
        map_int(modules, length) == 3
    )
  
  module_info <- 
    milnet %>%
    unnest(modules) %>% 
    transmute(
      module_id = modules,
      a0 = ifelse(module_id == "A1" | (double_repr & !grepl("^.1$", module_id)), 1, 0),
      burn = a0 > 0
    )
  
  module_network <- 
    map_df(seq_len(nrow(milnet)), function(i) {
      from_ <- milnet$from[[i]]
      mods <- milnet$modules[[i]]
      net <- 
        tibble(
          from = mods[-length(mods)],
          to = mods[-1]
        )
      if (from_ != "sBurn") {
        net <- bind_rows(
          tibble(
            from = milnet %>% filter(to == from_) %>% pull(modules) %>% first() %>% last(),
            to = mods[[1]]
          ),
          net
        )
      }
      net
    }) %>% 
    left_join(
      module_info %>% transmute(
        to = module_id,
        effect = ifelse(a0 > 0, -1, 1),
        strength = ifelse(a0 > 0, 10, 1),
        cooperativity = 2
      ), by = "to")
  
  module_network <- 
    inner_join(
      milnet %>% transmute(buff = from, from = map_chr(modules, first)),
      milnet %>% transmute(buff = from, to = map_chr(modules, first)),
      by = "buff"
    ) %>% 
    group_by(buff) %>% 
    filter(n() > 1) %>% 
    ungroup() %>% 
    select(-buff) %>% 
    mutate(
      effect = ifelse(from != to, -1, 1),
      strength = ifelse(from != to, 20, 1),
      cooperativity = 2
    ) %>% 
    bind_rows(module_network)
  
  expression_patterns <- 
    milnet %>% 
    transmute(
      from,
      to,
      module_progression = map_chr(modules, ~ paste0("+", ., collapse = ",")),
      start = from == "sBurn",
      burn = from == "sBurn",
      time = map_int(modules, length)
    )
  
  # add burn modules (diplicates are okay)
  expression_patterns$module_progression[[1]] <- 
    paste0(
      expression_patterns$module_progression[[1]],
      paste0(",+", module_info %>% filter(burn) %>% pull(module_id), collapse = "")
    )
  
  backbone(module_info, module_network, expression_patterns)
}

#' @export
#' @rdname backbone
backbone_binary_tree <- function(
  num_modifications = rbinom(1, size = 6, 0.25) + 1
) {
  backbone_branching(num_modifications = num_modifications, min_degree = 3, max_degree = 3)
}


#' @export
#' @rdname backbone
backbone_consecutive_bifurcating <- function() {
  backbone_branching(num_modifications = 2, max_degree = 3)
}

#' @export
#' @rdname backbone
backbone_trifurcating <- function() {
  backbone_branching(num_modifications = 1, min_degree = 4, max_degree = 4)
}

#' @export
#' @rdname backbone
backbone_converging <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
    "A1", 1, TRUE,
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
    "E3", 0, FALSE
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "A1", "B1", 1, 1, 2,
    "A1", "C1", 1, 1, 2,
    "B1", "C1", -1, 100, 2,
    "C1", "B1", -1, 100, 2, 
    "B1", "B2", 1, 1, 2,
    "B2", "B3", 1, 1, 2,
    "C1", "C2", 1, 1, 2,
    "C2", "C3", 1, 1, 2,
    "B3", "D1", 1, 1, 2,
    "C3", "D1", 1, 1, 2,
    "D1", "B1", -1, 100, 2,
    "D1", "C1", -1, 100, 2, 
    "D1", "B3", 1, 1, 2,
    "D1", "C3", 1, 1, 2,
    "D1", "D2", -1, 10, 2,
    "D2", "E1", -1, 1, 2,
    "B2", "E1", -1, 1, 2,
    "C2", "E1", -1, 1, 2,
    "E1", "E2", 1, 1, 2,
    "E2", "E3", 1, 1, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "sBurn1", "sBurn2", "+A1,+D2", TRUE, TRUE, 2,
    "sBurn2", "preB", "+B1,+B2", FALSE, TRUE, 2,
    "sBurn2", "preC", "+C1,+C2", FALSE, TRUE, 2,
    "preB", "sB", "+B3", FALSE, FALSE, 2,
    "preC", "sC", "+C3", FALSE, FALSE, 2,
    "sB", "sD", "+D1,+C1,+C2,+C3", FALSE, FALSE, 3,
    "sC", "sD", "+D1,+B1,+B2,+B3", FALSE, FALSE, 3,
    "sD", "sE", "+E1,+E2,+E3", FALSE, FALSE, 3
  )
  
  backbone(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname backbone
backbone_cycle <- function() {
  module_info <- tribble(
    ~module_id, ~a0, ~burn,
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
#' @rdname backbone
backbone_linear <- function(num_modules = 5) {
  module_info <- bind_rows(
    tibble(module_id = "A1", a0 = 1, burn = TRUE),
    tibble(module_id = paste0("B", seq_len(num_modules-1)), a0 = 0, burn = FALSE)
  )
  
  module_network <- tibble(
    from = module_info$module_id[-num_modules],
    to = module_info$module_id[-1],
    effect = 1,
    strength = 1,
    cooperativity = 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "sBurn", "sA", "+A1", TRUE, TRUE, 2,
    "sA", "sB", paste(paste0("+B", seq_len(num_modules-1)), collapse = ","), FALSE, FALSE, num_modules - 1
  )
  
  backbone(module_info, module_network, expression_patterns)
}


#' @export
#' @rdname backbone
backbone_disconnected <- function(
  
) {
  backbones <- list_backbones()
  backbones <- backbones[names(backbones) != "disconnected"]
  
  left_name <- sample(names(backbones), 1)[[1]]
  right_name <- sample(names(backbones), 1)[[1]]
  
  left <- backbones[[left_name]]()
  right <- backbones[[right_name]]()
  
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
  lmis <- lmi %>% filter(a0 > 0) %>% pull(module_id)
  rmis <- rmi %>% filter(a0 > 0) %>% pull(module_id)
  leps <- lep %>% filter(start) %>% pull(from)
  reps <- rep %>% filter(start) %>% pull(from)
  
  assert_that(
    length(lmis) > 0,
    length(rmis) > 0,
    length(leps) == 1,
    length(reps) == 1
  )
  
  common_modules <- paste0("D", 1:5)
  
  module_info <- bind_rows(
    tribble(
      ~module_id, ~a0, ~burn,
      "A1", 1, TRUE,
      "B1", 0, TRUE,
      "C1", 0, TRUE
    ),
    lmi %>% mutate(a0 = 0),
    rmi %>% mutate(a0 = 0),
    tibble(module_id = common_modules, a0 = 0, burn = FALSE)
  ) %>% select(-color)
  
  module_network <- bind_rows(
    tribble(
      ~from, ~to, ~effect, ~strength, ~cooperativity,
      "A1", "B1", 1, 1, 2,
      "A1", "C1", 1, 1, 2,
      "B1", "B1", 1, 1, 2,
      "B1", "C1", -1, 100, 2,
      "C1", "C1", 1, 1, 2,
      "C1", "B1", -1, 100, 2
    ),
    lmi %>% filter(a0 > 0) %>% transmute(from = "B1", to = module_id, effect = 1, strength = a0, cooperativity = 2),
    rmi %>% filter(a0 > 0) %>% transmute(from = "C1", to = module_id, effect = 1, strength = a0, cooperativity = 2),
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
        "sBurn", "sA", paste0("+", c("A1", common_modules), collapse = ","), TRUE, TRUE, 2,
        "sA", leps, "+B1", FALSE, TRUE, 2,
        "sA", reps, "+C1", FALSE, TRUE, 2
      ),
      lep %>% mutate(start = FALSE),
      rep %>% mutate(start = FALSE)
    )
  
  out <- backbone(module_info, module_network, expression_patterns)
  out$left_backbone <- left_name
  out$right_backbone <- right_name
  out
}