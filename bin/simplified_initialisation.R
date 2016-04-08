library(dyngen)

# formule testje
g <- 21
regs <- c(3,4,5)
inputs <- fsum(lapply(regs, function(r) fvar("x", r) / fvar("kg", r, g)))
numerator <- fvar("a0g", g) + (fvar("a1") * inputs)
denominator <- fcon(1) + inputs
production.form <- fvar("r", g) * numerator / denominator
production.form


# barabasi albert testje
amnt.genes <- 50
amnt.edges <- 200
ba.network <- generate.ba(amnt.nodes = amnt.genes, amnt.edges = amnt.edges, offset.exponent = 1.5)

formulae <- unlist(recursive = F, lapply(seq_len(amnt.genes), function(g) {
  
  # generate production formula and nu
  if (length(regs) == 0) {
    prod.formula <- fvar("r", g) * fvar("a0g", g)
  } else {
    inputs <- fsum(lapply(regs, function(r) fvar("x", r) / fvar("kg", r, g)))
    numerator <- fvar("a0g", g) + (fvar("a1") * inputs)
    denominator <- fcon(1) + inputs
    prod.formula <- fvar("r", g) * numerator / denominator
  }
  
  prod.nu <- rep(0, amnt.genes)
  prod.nu[[g]] <- 1
  
  # generate decay formula and nu
  decay.formula <- fvar("d", g) * fvar("x", g)
  
  decay.nu <- rep(0, amnt.genes)
  decay.nu[[g]] <- -1
  
  # return formulae
  list(
    list(formula = prod.formula, nu = prod.nu),
    list(formula = decay.formula, nu = decay.nu)
  )
}))

# formulae <- unlist(recursive = F, lapply(seq_len(amnt.genes), function(g) {
#   list(
#     generate.production.formula(g, ba.network$incoming.edges[[g]], amnt.genes),
#     generate.decay.formula(g, amnt.genes)
#   )
# }))

formulae.strings <- sapply(formulae, function(fl) fl$formula@string)
formulae.nus <- sapply(formulae, function(fl) fl$nu)

# 
# 
# ## generating the (initial) parameters of the system
# G <- seq_len(amnt.genes)
# R = setNames(rlnorm(length(G), log(10)), paste0("r", G))
# D = setNames(rep(1,length(G)), paste0("d", G))
# 
# K = sapply(kterms, function(kterm) R[paste0("r",str_replace(kterm, "k(\\d*)g\\d*", "\\1"))]/2)
# names(K) = kterms
# A0 = rep(0.1, length(G))
# names(A0) = paste0("a0g", G)
# 
# nu <- cbind2(diag(nrow=length(G), ncol=length(G)), -diag(nrow=length(G), ncol=length(G)))
# X0 = round(R/2, 0)
# names(X0) = paste0("x", G)
# 
