library(tidyverse)
devtools::load_all(".")

set.seed(1)
model <- 
  initialise_model(
    modulenet = modulenet_bifurcating_converging(),
    platform = platform_simple(n_cells = 1000, n_features = 5 * 10, pct_main_features = 1),
    tfgen_params = tfgen_random(percentage_tfs = 1, min_tfs_per_module = 3),
    simulation_params = simulation_default(total_time = 10, num_simulations = 32),
    simulation_setup = simulation_setup_custom(),
    verbose = TRUE,
    num_cores = 8
  ) %>% 
  generate_tf_network() %>% 
  generate_feature_network() %>% 
  generate_simulation_setup() %>%
  simulate_cells()

plot_module_network(model)
plot_feature_network(model)
plot_simulations(model)

model <- model %>% generate_goldstandard()

plot_gold_simulations(model)







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
