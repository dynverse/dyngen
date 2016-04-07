library(igraph)
library(dyngen)

# formule testje
g <- 21
regs <- c(3,4,5)
inputs <- fsum(lapply(regs, function(r) fvar("x", r) / fvar("kg", r, g)))
numerator <- fvar("a0g", g) + (fvar("a1") * inputs)
denominator <- 1 + inputs
production.form <- fvar("r", g) * numerator / denominator
production.form


# barabasi albert testje
amnt.genes <- 50
amnt.edges <- amnt.genes*3

# generate normal BA network, no modularity
net <- generate.ba(amnt.nodes = amnt.genes, amnt.edges = amnt.edges, reverse.edges = T, offset.exponent = 1.5)
# plot.igraph(graph_from_data_frame(net$data.frame), layout = layout_nicely, vertex.color = "white", vertex.size = 8, vertex.label.cex = 1)

# generate BA network with modularity
net <- generate.ba.with.modules(amnt.nodes = amnt.genes, amnt.edges = amnt.edges, reverse.edges = T, exp1 = .5, exp2 = 6)
plot.igraph(graph_from_data_frame(net$data.frame), layout = layout_nicely, vertex.color = "white", vertex.size = 8, vertex.label.cex = 1)
# plot.igraph(graph_from_data_frame(net$data.frame), layout = layout_with_kk, vertex.color = "white", vertex.size = 8, vertex.label.cex = 1)
plot(density(net$out.degree))

qplot(rank(-net$degree), net$degree) + theme_classic() + scale_y_log10()
# pheatmap(sapply(seq_len(amnt.genes), function(i) sapply(seq_len(amnt.genes), function(j) .jacc(net$incoming.edges[[i]], net$incoming.edges[[j]])))) # jaccard of incoming edges
# pheatmap(sapply(seq_len(amnt.genes), function(i) sapply(seq_len(amnt.genes), function(j) .jacc(net$outgoing.edges[[i]], net$outgoing.edges[[j]])))) # jaccard of outgoing edges

formulae <- unlist(recursive = F, lapply(seq_len(amnt.genes), function(g) {
  regs <- net$incoming.edges[[g]]
  
  # generate production formula and nu
  if (length(regs) == 0) {
    prod.formula <- fvar("r", g) * fvar("a0g", g)
  } else {
    inputs <- fsum(lapply(regs, function(r) fvar("x", r) / fvar("kg", r, g)))
    numerator <- fvar("a0g", g) + (fvar("a1") * inputs)
    denominator <- 1 + inputs
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

formulae.strings <- sapply(formulae, function(fl) fl$formula@string)
formulae.nus <- sapply(formulae, function(fl) fl$nu)
formulae.vars <- unique(unlist(lapply(formulae, function(fl) extract.variables(fl$formula))))
formulae.var.types <- sapply(formulae.vars, function(va) va@name)
xl <- tapply(formulae.vars, formulae.var.types, identity, simplify = F)

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
