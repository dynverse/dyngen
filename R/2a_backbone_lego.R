backbone_lego <- function(...) {
  mods <- list(...)
  
  backbone(
    module_info = map_df(mods, "module_info"),
    module_network = map_df(mods, "module_network"),
    expression_patterns = map_df(mods, "expression_patterns")
  )
}

backbone_lego_linear <- function(
  from_name, 
  to_name, 
  type = sample(c("simple", "double_repression", "flipflop"), 1),
  num_modules = c(simple = 2, double_repression = 3, flipflop = 4)[[type]]
) {
  assert_that(num_modules >= 1)
  
  module_ids <- c(
    paste0(from_name, seq_len(num_modules)),
    paste0(to_name, 1)
  )
  
  module_info <- tibble(
    module_id = module_ids %>% head(-1),
    a0 = 0,
    burn = FALSE
  )
  
  if (type == "simple") {
    module_network <- 
      tibble(
        from = module_ids %>% head(-1),
        to = module_ids %>% tail(-1),
        effect = 1,
        strength = 1,
        cooperativity = 2
      )
  } else if (type == "double_repression") {
    assert_that(num_modules >= 2)
    
    a0s <- c(0, rep(1, num_modules-1), 0)
    if (num_modules %% 2 == 1) {
      a0s[[1]] <- 0
    }
    module_info <- module_info %>% mutate(
      a0 = a0s %>% head(-1),
      burn = a0 > 0
    )
    
    module_network <- 
      tibble(
        from = module_ids %>% head(-1),
        to = module_ids %>% tail(-1),
        effect = ifelse(a0s[-1] == 1, -1, 1),
        strength = ifelse(effect == 1, 1, 10),
        cooperativity = 2
      )
  } else if (type == "flipflop") {
    assert_that(num_modules >= 4)
    
    module_info <- tibble(
      module_id = module_ids %>% head(-1),
      a0 = 0,
      burn = FALSE
    )
    
    module_network <- bind_rows(
      tibble(
        from = module_ids %>% head(-1),
        to = module_ids %>% tail(-1),
        effect = 1,
        strength = 1,
        cooperativity = 2
      ),
      tribble(
        ~from, ~to, ~effect, ~strength, ~cooperativity,
        nth(module_ids, -2), nth(module_ids, 1), -1, 10, 2,
        nth(module_ids, -3), nth(module_ids, -1), -1, 100, 2,
        nth(module_ids, -2), nth(module_ids, -2), 1, 1, 2
      )
    )
  }
  
  expression_patterns <- tibble(
    from = paste0("s", from_name),
    to = paste0("s", to_name),
    module_progression = paste0("+", module_ids %>% head(-1), collapse = ","),
    start = FALSE,
    burn = FALSE,
    time = case_when(
      type == "flipflop" ~ num_modules * 2,
      TRUE ~ num_modules
    )
  )
  
  lst(module_info, module_network, expression_patterns)
}

backbone_lego_branching <- function(
  from_name, 
  to_name, 
  type = "simple",
  num_modules = length(to_name) + 1
) {
  assert_that(
    length(to_name) >= 2,
    num_modules >= length(to_name) + 1
  )
  
  lin_length <- num_modules - length(to_name)
  
  my_module_ids <- 
    paste0(from_name, seq_len(lin_length))
  our_module_ids <- 
    paste0(from_name, seq(lin_length + 1, num_modules))
  their_module_ids <- 
    paste0(to_name, 1)
  module_ids <- c(my_module_ids, our_module_ids, their_module_ids)
  
  module_info <- tibble(
    module_id = c(my_module_ids, our_module_ids),
    a0 = 0,
    burn = FALSE
  )
  
  module_network <- bind_rows(
    tibble(
      from = my_module_ids %>% head(-1),
      to = my_module_ids %>% tail(-1),
      effect = 1,
      strength = 1,
      cooperativity = 2,
    ),
    tibble(
      from = last(my_module_ids),
      to = our_module_ids,
      effect = 1,
      strength = 1,
      cooperativity = 2
    ),
    tibble(
      from = our_module_ids,
      to = our_module_ids,
      effect = 1,
      strength = 5,
      cooperativity = 2
    ),
    tibble(
      from = our_module_ids,
      to = first(my_module_ids),
      effect = -1,
      strength = 5,
      cooperativity = 2
    ),
    tibble(
      from = our_module_ids,
      to = their_module_ids,
      effect = 1,
      strength = 1,
      cooperativity = 2
    ),
    crossing(
      from = our_module_ids,
      to = our_module_ids
    ) %>% 
      filter(from != to) %>% 
      mutate(
        effect = -1,
        strength = 100,
        cooperativity = 2
      )
  )
  
  in_edge <- if (length(my_module_ids) >= 1) {
    paste0("+", my_module_ids, ",", collapse = "")
  } else {
    c()
  }
  
  expression_patterns <- tibble(
    from = paste0("s", from_name),
    to = paste0("s", to_name),
    module_progression = paste0(in_edge, "+", our_module_ids),
    start = FALSE,
    burn = FALSE,
    time = num_modules + 1
  )
  
  lst(module_info, module_network, expression_patterns)
}


backbone_lego_burn <- dynutils::inherit_default_params(
  list(backbone_lego_linear), 
  function(to_name, type, num_modules) {
    out <- backbone_lego_linear("Burn", to_name, type = type, num_modules = num_modules)
    out$module_info <- out$module_info %>% mutate(
      a0 = ifelse(module_id == "Burn1", 1, a0),
      burn = TRUE
    )
    out$expression_patterns <- out$expression_patterns %>% mutate(
      start = TRUE,
      burn = TRUE
    )
    out
  }
)

backbone_lego_end <- dynutils::inherit_default_params(
  list(backbone_lego_linear), 
  function(from_name, type, num_modules) {
    out <- backbone_lego_linear(from_name, paste0("End", from_name), type = type, num_modules = num_modules)
    out$module_network <- out$module_network %>% filter(!grepl("^End", to))
    out$expression_patterns <- out$expression_patterns %>% mutate(
      time = time - 1
    ) %>% 
      filter(time > 0)
    out
  }
)