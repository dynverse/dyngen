generate_kinetics = function(vargroups, variables, nus.changes) {
  ## generating the (initial) parameters of the system
  
  R = sapply(vargroups$r, function(r) 20)
  D = sapply(vargroups$d, function(d) 5)
  
  K = sapply(vargroups$k, function(k) unname((R/D)[paste0("r_",str_replace(k, "k_\\d*_(\\d*)", "\\1"))]/2/variables[[k]]$strength))
  
  P = sapply(vargroups$p, function(p) 0.2)
  Q = sapply(vargroups$q, function(q) 0.2)
  
  g = c(1)
  #P[paste0("p_", g)] = P[paste0("p_", g)] /5
  #Q[paste0("q_", g)] = Q[paste0("q_", g)] /5
  
  A0 = map_dbl(variables[vargroups$a0], "value")
  #A0[[1]] = 1
  #A0[[4]] = 0
  #A0[paste0("a0_", ldtfs)] = 1 # IF CYCLE
  #A0[modules[[1]]] = 1
  
  A = sapply(vargroups$a, function(a) ifelse(variables[[a]]$effect==1, 1, 0))
  
  initial.state = numeric()
  initial.state[vargroups$x] = 0
  initial.state[vargroups$y] = 0
  initial.state[vargroups$u] = 1
  initial.state[vargroups$b] = 0
  #initial.state["rg"] = 1
  molecules = names(initial.state)
  
  formulae.nus = sapply(nus.changes, function(nu.changes) {
    nu = sapply(molecules, function(r) 0)
    nu[names(nu.changes)] = nu.changes
    nu
  })
  
  print(molecules)
  
  params = c(R, D, K, P, Q, A0,A, rg=1,dg=1,kg=1,pg=10,qg=10, a1=1)
  
  named.list(formulae.nus, params, initial.state)
}