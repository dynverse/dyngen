#' Design your own custom backbone easily
#' 
#' You can use the `bblego` functions in order to create 
#' custom backbones using various components. Please note that the `bblego` 
#' functions currently only allow you to create tree-like backbones. 
#' 
#' A backbone always needs to start with a single `bblego_start()` state and
#' needs to end with one or more `bblego_end()` states.
#' The order of the mentioned states needs to be such that a state is never
#' specified in the first argument (except for `bblego_start()`) before
#' having been specified as the second argument.
#' 
#' @param from The begin state of this component.
#' @param to The end state of this component.
#' @param type Some components have alternative module regulatory networks. 
#' 
#' `bblego_start()`, `bblego_linear()`, `bblego_end()`: 
#' 
#' * `"simple"`: a sequence of modules in which every module upregulates the next module.
#' * `"doublerep1"`: a sequence of modules in which every module downregulates the next
#'   module, and each module has positive basal expression.
#' * `"doublerep2"`: a sequence of modules in which every module upregulates the next 
#'   module, but downregulates the one after that.
#' * `"flipflop"`: a sequence of modules in which every module upregulates the next module.
#'   In addition, the last module upregulates itself and strongly downregulates the first module.
#'   
#' `bblego_branching()`: 
#' 
#' * `"simple"`: a set of `n` modules (with `n = length(to)`) which all downregulate
#'   one another and upregulate themselves. This causes a branching to occur in the trajectory.
#'   
#' @param num_modules The number of modules this component is allowed to use. 
#'   Various components might require a minimum number of components in order to work properly.
#' @param burn Whether or not these components are part of the warm-up simulation.
#' @param ...,.list `bblego` components, either as separate args or as a list.
#' 
#' @export
#' 
#' @examples 
#' backbone <- bblego(
#'   bblego_start("A", type = "simple", num_modules = 2),
#'   bblego_linear("A", "B", type = "simple", num_modules = 3),
#'   bblego_branching("B", c("C", "D"), type = "simple", num_steps = 3),
#'   bblego_end("C", type = "flipflop", num_modules = 4),
#'   bblego_end("D", type = "doublerep1", num_modules = 7)
#' )
bblego <- function(..., .list = NULL) {
  mods <- c(list(...), .list)
  
  backbone(
    module_info = map_df(mods, "module_info"),
    module_network = map_df(mods, "module_network"),
    expression_patterns = map_df(mods, "expression_patterns")
  )
}

#' @rdname bblego
#' @export
bblego_linear <- function(
  from, 
  to, 
  type = sample(c("simple", "doublerep1", "doublerep2"), 1),
  num_modules = sample(4:6, 1),
  burn = FALSE
) {
  assert_that(num_modules >= 1)
  
  module_ids <- c(
    paste0(from, seq_len(num_modules)),
    paste0(to, 1)
  )
  
  module_info <- tibble(
    module_id = module_ids %>% head(-1),
    basal = 0,
    burn = burn,
    independence = 1
  )
  
  if (type == "simple") {
    module_network <- 
      tibble(
        from = module_ids %>% head(-1),
        to = module_ids %>% tail(-1),
        effect = 1L,
        strength = 1,
        hill = 2
      )
  } else if (type == "doublerep1") {
    assert_that(num_modules >= 2)
    
    bas <- c(0, rep(1, num_modules-1), 0)
    if (num_modules %% 2 == 0) {
      bas[[2]] <- 0
    }
    module_info <- module_info %>% mutate(
      basal = bas %>% head(-1),
      burn = burn | basal > 0
    )
    
    module_network <- 
      tibble(
        from = module_ids %>% head(-1),
        to = module_ids %>% tail(-1),
        effect = ifelse(bas[-1] > 0, -1L, 1L),
        strength = ifelse(effect == 1, 1, 10),
        hill = 2
      )
  } else if (type == "doublerep2") {
    assert_that(num_modules >= 4)
    pos_to <- 
      if (num_modules %% 2 == 1) {
        module_ids[c(2, num_modules, num_modules + 1)]
      } else {
        module_ids[c(2, num_modules + 1)]
      }
    module_network <- 
      bind_rows(
        tibble(
          from = module_ids %>% head(-1),
          to = module_ids %>% tail(-1),
          effect = ifelse(to %in% pos_to, 1L, -1L),
          strength = 1,
          hill = 2
        ),
        tibble(
          from = module_ids %>% head(-2) %>% head(-1),
          to = module_ids %>% tail(-2) %>% head(-1),
          effect = 1L,
          strength = seq_along(from),
          hill = 2
        ) %>% 
          head(., ifelse(num_modules %% 2 == 1, nrow(.) - 1, nrow(.)))
      )
  } else if (type == "flipflop") {
    assert_that(num_modules >= 4)
    
    module_info <- tibble(
      module_id = module_ids %>% head(-1),
      basal = 0,
      burn = burn,
      independence = 1
    )
    
    module_network <- bind_rows(
      tibble(
        from = module_ids %>% head(-1),
        to = module_ids %>% tail(-1),
        effect = 1L,
        strength = 1,
        hill = 2
      ),
      tribble(
        ~from, ~to, ~effect, ~strength, ~hill,
        nth(module_ids, -2), nth(module_ids, 1), -1L, 10, 2,
        nth(module_ids, -3), nth(module_ids, -1), -1L, 100, 2,
        nth(module_ids, -2), nth(module_ids, -2), 1L, 1, 2
      )
    )
  } else {
    stop("Unknown 'type' parameter")
  }
  
  expression_patterns <- tibble(
    from = paste0("s", from),
    to = paste0("s", to),
    module_progression = paste0("+", module_ids %>% head(-1), collapse = ","),
    start = FALSE,
    burn = burn,
    time = case_when(
      type == "flipflop" ~ num_modules * 60,
      num_modules == 1 ~ 40,
      num_modules == 2 ~ 60,
      TRUE ~ as.numeric(num_modules) * 30
    )
  )
  
  lst(module_info, module_network, expression_patterns)
}



#' @rdname bblego
#' @param num_steps The number of branching steps to reduce the odds of double positive cells occurring. 
#' @export
bblego_branching <- function(
  from,
  to,
  type = "simple",
  num_steps = 3,
  num_modules = 2 + length(to) * (3 + num_steps),
  burn = FALSE
) {
  reversible <- FALSE
  
  assert_that(
    length(to) >= 2,
    num_modules >= (3 + num_steps) * length(to) + 2
  )

  type <- match.arg(type)

  num_start <- num_modules - (3 + num_steps) * length(to) - 1
  starts <- paste0("START", seq_len(num_start)-1)
  xor <- "XOR"
  me <- paste0("step", to, "0") %>% set_names(to)
  notme <- paste0("not", to) %>% set_names(to)
  meandnotrest <- paste0(to, "andnotrest") %>% set_names(to)
  step_groups <- map(to, ~ paste0("step", ., seq_len(num_steps))) %>% set_names(to)
  their <- paste0(to, "1")

  module_info <- tibble(
    tmp = c(starts, me, notme, meandnotrest, xor, unlist(step_groups)),
    module_id = paste0(from, seq_along(tmp)),
    basal = ifelse(tmp %in% notme, 1, 0),
    burn = basal > 0,
    independence = ifelse(tmp %in% c(meandnotrest, unlist(step_groups)), 0, ifelse(tmp %in% me, .5, 1))
  )
  name_map <- c(module_info %>% select(tmp, module_id) %>% deframe(), set_names(their, their))
  module_info <- module_info %>% select(-tmp)

  module_network <- bind_rows(
    modnet_chain(starts),
    modnet_edge(last(starts), me),
    modnet_pairwise(me, effect = -1L, strength = 10),
    
    modnet_edge(me, notme, effect = -1L, strength = 100),
    modnet_edge(me, meandnotrest),
    modnet_pairwise(notme, meandnotrest),
    modnet_edge(meandnotrest, xor),
    map_df(to, function(too) {
      bind_rows(
        modnet_edge(me[[too]], step_groups[[too]]),
        modnet_edge(xor, step_groups[[too]])
      )
    }),
    modnet_edge(map_chr(step_groups, last), their)
  )
  if (!reversible) {
    module_network <- module_network %>% bind_rows(
      modnet_self(me, effect = 1L, strength = 2),
    )
  }
  
  module_network <- module_network %>% 
    mutate(
      from = name_map[from], 
      to = name_map[to],
      hill = 4 # override hill parameter to counter the repeated AND gates
    )

  expression_patterns <- bind_rows(
    tibble(
      from = paste0("s", from),
      to = paste0("s", !!from, "mid"),
      module_progression = paste0("+", name_map[c(starts, notme)], collapse = ","),
      start = FALSE,
      burn = burn,
      time = (length(starts) + 1) * 20
    ),
    tibble(
      from = paste0("s", from, "mid"),
      to = paste0("s", to),
      module_progression = map_chr(!!to, ~ paste0("+", name_map[c(me[[.]], meandnotrest[[.]], xor, step_groups[[.]])], collapse = ",")),
      start = FALSE,
      burn = burn,
      time = (num_steps + 1) * 30
    )
  )

  lst(module_info, module_network, expression_patterns)
}

#' @rdname bblego
#' @export
bblego_start <- dynutils::inherit_default_params(
  list(bblego_linear), 
  function(to, type, num_modules) {
    out <- bblego_linear(
      from = "Burn", 
      to = to,
      type = type,
      num_modules = num_modules,
      burn = TRUE
    )
    out$module_info <- 
      out$module_info %>% 
      mutate(basal = ifelse(module_id == "Burn1", 1, basal))
    out$expression_patterns <-
      out$expression_patterns %>% 
      mutate(start = TRUE)
    out
  }
)

#' @rdname bblego
#' @export
bblego_end <- dynutils::inherit_default_params(
  list(bblego_linear), 
  function(from, type, num_modules) {
    out <- bblego_linear(from, paste0("End", from), type = type, num_modules = num_modules)
    out$module_network <- out$module_network %>% filter(!grepl("^End", to))
    out$expression_patterns <- out$expression_patterns %>% mutate(
      time = time - 1
    ) %>% 
      filter(time > 0)
    out
  }
)