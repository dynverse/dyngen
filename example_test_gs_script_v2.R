library(tidyverse)
devtools::load_all(".")

model <- read_rds("~/yay3.rds")

plot_module_network(model)

sim_params <- model$simulation_params
sim_system <- model$simulation_system

propensity_funs <-
  sim_system$formulae %>% 
  select(formula_id, formula) %>% 
  deframe


# find gold transition order
mod_changes <- 
  model$modulenet$expression_patterns %>% 
  mutate(
    mod_diff = module_progression %>% str_split("\\|"),
    substate = map(mod_diff, seq_along)
  ) %>% 
  select(-module_progression) %>% 
  unnest(mod_diff, substate) %>% 
  mutate(
    mod_diff = str_split(mod_diff, ","),
    mod_on = map(mod_diff, function(x) x %>% keep(grepl("\\+", .)) %>% gsub("\\+", "", .)),
    mod_off = map(mod_diff, function(x) x %>% keep(grepl("-", .)) %>% gsub("-", "", .))
  ) %>% 
  group_by(from, to) %>% 
  mutate(
    from_ = ifelse(row_number() == 1, from, c("", paste0(from, to, letters[row_number() - 1]))),
    to_ = ifelse(row_number() == n(), to, from_[row_number() + 1])
  ) %>% 
  ungroup()

gr <- igraph::graph_from_data_frame(mod_changes %>% select(from_, to_))

mod_order <- igraph::topo_sort(gr) %>% names()


# start constructing golds
gold_sim_outputs <- list()
gold_sim_vectors <- list()
gold_sim_modson <- list()

gold_sim_vectors[[mod_order[[1]]]] <- sim_system$initial_state %>% as.matrix
gold_sim_modson[[mod_order[[1]]]] <- c()

for (i in seq(2, length(mod_order))) {
  to_ <- mod_order[[i]]
  
  for (j in which(mod_changes$to_ == to_)) {
    from_ <- mod_changes$from_[[j]]
    mod_on <- mod_changes$mod_on[[j]]
    mod_off <- mod_changes$mod_off[[j]]
    prev_mods <- gold_sim_modson[[from_]]
    mods <- union(mod_on, prev_mods) %>% setdiff(mod_off)
    
    cat("Processing ", to_, " := ", from_, " ~ ", paste(mods, collapse = "|"), "\n", sep = "")
    
    # fetch initial state
    new_initial_state <- rowMeans(gold_sim_vectors[[from_]])

    # generate nus
    new_nus <- sim_system$nus
    filt <- colnames(new_nus) %>% purrr:::discard(grepl(paste0("_(", paste(mods, collapse = "|"), ")_"), .))
    new_nus[, filt] <- new_nus[, filt] * 0

    # actual simulation
    out <- fastgssa::ssa(
      initial.state = new_initial_state,
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
    simulation <-
      .simulate_cells_process_ssa(out) %>%
      mutate(
        j,
        from_,
        to_
      ) %>%
      select(j, from_, to_, everything())

    end_state <-
      simulation %>%
      slice(n()) %>%
      select(one_of(sim_system$molecule_ids)) %>%
      .[1, , drop = TRUE] %>%
      unlist() %>%
      as.matrix()

    gold_sim_outputs[[j]] <- simulation
    
    gold_sim_modson[[to_]] <- mods
    
    if (!to_ %in% gold_sim_vectors) {
      gold_sim_vectors[[to_]] <- end_state
    } else {
      gold_sim_vectors[[to_]] <- cbind(gold_sim_vectors[[to_]], end_state)
    }
    
  }
}
  

########################################
###             combine              ###
########################################

simulation <- 
  bind_rows(gold_sim_outputs) %>% 
  left_join(mod_changes %>% select(from, to, from_, to_, substate), by = c("from_", "to_")) %>% 
  group_by(from, to) %>% 
  mutate(
    simulation_i = 0,
    type = paste0(from, "->", to),
    time = ((substate - 1) * 2 + t) / 2 / max(substate)
  ) %>% 
  ungroup()

model_new <- model
model_new$simulations <- bind_rows(simulation, model$simulations %>% mutate(type = paste0("simulation")))

########################################
###              plot                ###
########################################
# expr <- model_new$simulations %>% select(one_of(model_new$simulation_system$molecule_ids)) %>% as.matrix
sim_f <- model_new$simulations %>% select(t, simulation_i, type, from_, to_, j, from, to, time, substate)
expr <- model_new$simulations %>% select(-one_of(colnames(sim_f))) %>% as.matrix

space <- SCORPIUS::reduce_dimensionality(expr, dist_fun = SCORPIUS::correlation_distance, num_landmarks = 2000, ndim = 10)
plot_df <- bind_cols(sim_f, as.data.frame(space))

g1 <- ggplot(mapping = aes(Comp1, Comp4, colour = type, group = paste0(type, "_", simulation_i))) +
  geom_path(data = plot_df %>% filter(simulation_i > 0), colour = "darkgray") +
  geom_path(data = plot_df %>% filter(simulation_i == 0), size = 2) +
  theme_bw()

g2 <- ggplot(mapping = aes(Comp1, Comp4, colour = paste0(from_, "->", to_), group = paste0(type, "_", simulation_i))) +
  geom_path(data = plot_df %>% filter(simulation_i > 0), colour = "darkgray") +
  geom_path(data = plot_df %>% filter(simulation_i == 0), size = 2) +
  theme_bw()

patchwork::wrap_plots(g1, g2)



ix <- sim_f$simulation_i == 0
sim_f_tr <- data.frame(sim_f[ix,], space[ix, ]) %>% mutate(group = "train")
expr_tr <- expr[ix,]

sim_f_pr <- data.frame(sim_f[-ix,], space[-ix, ]) %>% mutate(group = "pred")
expr_pr <- expr[-ix,]

rf <- ranger::ranger(type ~ ., data.frame(type = sim_f_tr$type, expr_tr, check.names = FALSE, stringsAsFactors = FALSE))
pred <- predict(rf, data.frame(expr_pr, check.names = FALSE, stringsAsFactors = FALSE))
sim_f_pr$pred <- pred$predictions


g <- ggplot(mapping = aes(Comp1, Comp4, group = paste0(type, "_", simulation_i))) +
  geom_path(aes(colour = pred), sim_f_pr, size = 2) +
  geom_path(aes(colour = type), sim_f_tr, size = 2) +
  theme_bw() +
  facet_wrap(~ simulation_i)
g
# ggsave("~/yay.pdf", g, width = 20, height = 16)

rf <- ranger::ranger(type ~ ., data.frame(type = sim_f_tr$type, expr_tr, check.names = FALSE, stringsAsFactors = FALSE), probability = TRUE)
pred <- predict(rf, data.frame(expr_pr, check.names = FALSE, stringsAsFactors = FALSE))$predictions

milestone_network <- 
  sim_f_tr %>% 
  group_by(type) %>%
  summarise(length = sqrt(sum(diff(Comp1)^2 + diff(Comp2)^2 + diff(Comp3)^3))) %>% 
  filter(type != "burn") %>% 
  transmute(
    from = gsub("->.*", "", type),
    to = gsub(".*->", "", type),
    length = length / sum(length) * 20,
    directed = TRUE
  )

rownames(expr_pr) <- rownames(pred) <- paste0("Cell_", seq_len(nrow(pred)))


rf <- randomForestSRC::rfsrc(Multivar(type, t) ~ ., data.frame(type = sim_f_tr$type, t = sim_f_tr$t, expr_tr, check.names = FALSE))
pred <- predict(rf, data.frame(expr_pr, check.names = FALSE, stringsAsFactors = FALSE))

progressions <- 
  tibble(
    cell_id = rownames(expr_pr),
    type = pred$classOutput$type$class,
    t = pred$regrOutput$t$predicted
  ) %>% 
  filter(type != "burn") %>% 
  transmute(
    cell_id,
    from = gsub("->.*", "", type),
    to = gsub(".*->", "", type),
    percentage = t / 2
  )

counts <- expr_pr[progressions$cell_id, model$feature_info$x] %>% round()
keep <- which(rowSums(counts) != 0) %>% names() %>% sample(2000)

progressions <- progressions %>% filter(cell_id %in% keep)
counts <- counts[keep, , drop = FALSE]
traj <- 
  wrap_expression(
    counts = counts,
    expression = counts
  ) %>% 
  add_trajectory(
    milestone_network = milestone_network,
    progressions = progressions
  )

g1 <- dynplot::plot_graph(traj)
dr <- dyndimred::dimred_landmark_mds(counts)
g2 <- dynplot::plot_dimred(traj, dimred = dr + rnorm(length(dr), 0, .01))
patchwork::wrap_plots(g1, g2)
