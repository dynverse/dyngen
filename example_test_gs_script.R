model <- read_rds("~/yay.rds")

plot_module_network(model)

sim_params <- model$simulation_params
sim_system <- model$simulation_system

set.seed(1)

propensity_funs <-
  sim_system$formulae %>% 
  select(formula_id, formula) %>% 
  deframe

nus_burn <- sim_system$nus
nus_burn[setdiff(rownames(nus_burn), sim_system$burn_variables),] <- 0

# burn in
out <- fastgssa::ssa(
  initial.state = sim_system$initial_state, 
  propensity.funs = propensity_funs,
  nu = nus_burn %>% as.matrix,
  final.time = sim_params$burn_time, 
  parms = sim_system$parameters,
  method = sim_params$ssa_algorithm,
  recalculate.all = TRUE, 
  stop.on.negstate = FALSE,
  stop.on.propensity = FALSE,
  verbose = FALSE
)

burn_out <- .simulate_cells_process_ssa(out) %>% mutate(type = "burn")

new_initial_state <- 
  burn_out %>% 
  slice(n()) %>% 
  select(one_of(sim_system$molecule_ids)) %>% 
  .[1, , drop = TRUE] %>% 
  unlist()


helper <- function(module, state, name) {
  new_nus <- sim_system$nus
  filt <- colnames(new_nus) %>% purrr:::discard(grepl(paste0("_", module, "_"), .))
  # new_nus[, filt] <- 0
  
  new_nus[, filt] <- new_nus[, filt] * .1
  
  # actual simulation
  out <- fastgssa::ssa(
    initial.state = state,
    propensity.funs = propensity_funs,
    nu = new_nus %>% as.matrix,
    final.time = sim_params$burn_time,
    parms = sim_system$parameters,
    method = sim_params$ssa_algorithm,
    recalculate.all = TRUE, 
    stop.on.negstate = FALSE,
    stop.on.propensity = FALSE,
    verbose = FALSE
  )
  
  # add both burnin as normal simulation together
  simulation <- .simulate_cells_process_ssa(out) %>% mutate(type = name)
  
  end_state <- 
    simulation %>% 
    slice(n()) %>% 
    select(one_of(sim_system$molecule_ids)) %>% 
    .[1, , drop = TRUE] %>% 
    unlist()
  
  lst(simulation, end_state)
}

########################################
###           transitions            ###
########################################

m2 <- helper(module = "M2", state = new_initial_state, "burn->m2")
m3 <- helper(module = "M3", state = m2$end_state, "m2->m3")
m4 <- helper(module = "M4", state = m2$end_state, "m2->m4")
m10 <- helper(module = "M10", state = m4$end_state, "m4->m10")
m5 <- helper(module = "M5", state = m3$end_state, "m3->m5")
m6 <- helper(module = "M6", state = m5$end_state, "m5->m6")
m7 <- helper(module = "M7", state = m5$end_state, "m5->m7")
m8 <- helper(module = "M8", state = m6$end_state, "m6->m8")
m9 <- helper(module = "M9", state = m7$end_state, "m7->m9")


########################################
###             combine              ###
########################################

simulation <- 
  bind_rows(
    burn_out,
    m2$simulation,
    m3$simulation,
    m4$simulation,
    m10$simulation,
    m5$simulation,
    m6$simulation,
    m7$simulation,
    m8$simulation,
    m9$simulation
  ) %>% 
  mutate(simulation_i = 0)

model_new <- model
model_new$simulations <- bind_rows(simulation, model$simulations %>% mutate(type = paste0("simulation")))

########################################
###              plot                ###
########################################
expr <- model_new$simulations %>% select(one_of(model_new$simulation_system$molecule_ids)) %>% as.matrix
expr <- expr[,colSums(expr) != 0]
sim_f <- model_new$simulations[rowSums(expr) != 0,]
expr <- expr[rowSums(expr) != 0,]
space <- SCORPIUS::reduce_dimensionality(expr, dist_fun = SCORPIUS::correlation_distance, num_landmarks = 2000)
plot_df <- bind_cols(sim_f %>% select(t, simulation_i, type), as.data.frame(space))

ggplot() +
  geom_path(aes(Comp1, Comp2, colour = type, group = paste0(type, "_", simulation_i)), plot_df %>% filter(simulation_i > 0), colour = "darkgray") +
  geom_path(aes(Comp1, Comp2, colour = type, group = paste0(type, "_", simulation_i)), plot_df %>% filter(simulation_i == 0), size = 2) +
  scale_color_brewer(palette = "Set3") +
  theme_bw()
