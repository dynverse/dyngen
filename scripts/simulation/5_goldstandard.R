get_module_counts = function(counts, modulemembership) {
  lapply(modulemembership, function(module) apply(counts[,module], 1, mean)) %>% do.call(cbind, .)
}


get_states_and_progression = function(expression_modules, expression_cutoff=3, expression_max=4) {
  expression_modules %>% apply(1, function(x) {
    state = (1:length(x))[x>expression_cutoff] %>% last
    if (is.na(state)) {
      warning("no good state")
      state = which.max(x)
    }
    progression = pmax(pmin((x[state] - expression_cutoff)/(expression_max-expression_cutoff), 1), 0)
    
    tibble(state=state, progression=progression)
  }) %>% bind_rows %>% mutate(cell=rownames(expression_modules))
}