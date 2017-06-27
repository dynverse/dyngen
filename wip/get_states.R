# shouldn't we start from the end and retrace our steps? What if due to stochsticity a primed cell still choses for another path?

get_states = function(expression_modules, model, minexpression = 0.8, minratio = 1.5, state=1) {
  statechanged = T
  states = state
  for (cell in rownames(expression_modules)[2:nrow(expression_modules)]) {
    if (statechanged) {
      nextstates = model$statenet %>% filter(from == state) %>% .$to
      state2module = model$modulenodes %>% filter(!is.na(state)) %>% {set_names(.$module, .$state)}
      nextmodules = state2module[nextstates %>% as.character]
    }
    
    statechanged = F
    if(length(nextmodules) == 1) {
      if(expression_modules[cell,nextmodules] > minexpression) {
        state = nextstates[[1]]
        statechanged = T
      }
    } else if (length(nextmodules) > 1) {
      expression_modules_cell = expression_modules[cell,nextmodules]
      if(any(expression_modules_cell > minexpression)) {
        ratios = expression_modules_cell[which.max(expression_modules_cell)]/expression_modules_cell[-which.max(expression_modules_cell)]
        if(all(ratios > minratio)) {
          state = nextstates[[which.max(expression_modules_cell)]]
          statechanged = T
        }
      }
      
    }
    states = c(states, state)
  }
  
  tibble(state=states, cell=rownames(expression_modules))
}

combined = matrices %>% do.call(rbind, .)
nrows = map_int(matrices, nrow) %>% c(0, .)
newmatrices = map(seq_len(length(nrows)-1), function(i) {combined[sum(nrows[1:i]):sum(nrows[1:(i+1)]), ]})
map(seq_len(length(nrows)), ~rep(., nrows[[.]])) %>% unlist

quant_scale_combined = function(matrices, simulations, outlier.cutoff = 0.05) {
  combined = map(simulations, "expression") %>% do.call(rbind2, .)
  combined_scaled = combined %>% SCORPIUS::quant.scale(outlier.cutoff=0.05)
  attributes(combined_scaled)
}
