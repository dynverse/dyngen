adapt_model = function(vars, model) {
  names(vars) = dynamic_params
  model$params = c(static_paramvalues, vars)
  
  model
}

run_cycle = function(model, totaltime=50, burntime=5) {
  experiment = dyngen:::simulate_multiple(model, burntime, totaltime,  ssa.algorithm=fastgssa::ssa.em(noise_strength=1), local=TRUE, nsimulations = 4)
  experiment$expression = map(experiment, "expression") %>% do.call(rbind, .)
  experiment$expression = experiment$expression[sample(seq_len(nrow(experiment$expression)), 1000),]
  experiment$cellinfo = map(experiment, "stepinfo") %>% bind_rows() %>% filter(stepid %in% rownames(experiment$expression)) %>% arrange(simulationtime)
  experiment$expression = experiment$expression[experiment$cellinfo$stepid, ]
  
  if(!is.null(model$modulemembership)) {experiment$expression_modules = dyngen:::get_module_counts(experiment$expression, model$modulemembership)}
  
  experiment
}

fitness_cycle = function(experiment) {
  plot(experiment$expression_modules[, "M1"])
  delayinfo = cal_delay(experiment$expression_modules[, "M1"])
  timedelay = experiment$cellinfo$simulationtime[delayinfo$peak]
  periodgoal = 10
  
  #print(timedelay)
  #print(delayinfo)
  
  -abs(periodgoal-timedelay) - (1-delayinfo$acf)
  
  #xs = experiment$cellinfo$simulationtime#seq(0, 30, by=0.1)
  #ys = sin((xs + pi/2)*pi*2/periodgoal)
  
  #plot(ys, experiment$expression_modules[, "M1"])
  #plot(c(xs, xs), c(ys*100, experiment$expression_modules[, "M1"]))
  
  #plot(smoother::smth(experiment$expression_modules[, "M1"], tails=T, window=200))
  #print(smoother::smth(experiment$expression_modules[, "M1"], tails=T, window=200) %>% find_peaks(100))
  
  #cor(ys, experiment$expression_modules[, "M1"])
}

fitness_cycle = function(experiment) {
  expression_smooth = experiment$expression_modules %>% zoo::rollmean(50, c("extend", "extend", "extend")) %>% set_rownames(rownames(experiment$expression))
  expression_smooth %>% apply(1, which.max) %>% rle() %>% .$values %>% paste0(collapse="") %>% str_count("123")
}

fitness_cycle = function(experiment) {
  expression_smooth = experiment$expression_modules %>% zoo::rollmean(50, c("extend", "extend", "extend")) %>% set_rownames(rownames(experiment$expression))
  
  delays = seq(1, 6, length.out=50)
  shifts = seq(0, 10, length.out=20)
  design = expand.grid(delay=delays, shift=shifts)
  cors = purrr::map_dbl(seq_len(nrow(design)), ~cor(sin((2*pi*experiment$cellinfo$simulationtime+design[., "shift"])/design[., "delay"]), expression_smooth[, "M1"]))
  #plot(cors)
  #print(design[which.max(cors),])
  
  cors %>% max
}

fitness_bifurcating = function(experiment) {
  modexpression = SCORPIUS::quant.scale(experiment$expression_modules)[, c("M2", "M3")]
  quantile_mean = modexpression %>% apply(2, function(x) quantile(x, 0.8)) %>% mean()
  diff_mean = modexpression %>% {.[, 1]-.[, 2]} %>% abs %>% mean
  
  quantile_mean * diff_mean
}

fit_cycle = function(params) {
  model = adapt_model(params, base_model)
  experiment = run_cycle(model)
  fitness_cycle(experiment)
}

fit_bifurcating = function(params) {
  model = adapt_model(params, base_model)
  experiment = run_cycle(model)
  fitness_cycle(experiment)
}