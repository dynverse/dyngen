library(fastgssa)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(pheatmap)
library(parallel)
library(PRISM)

library(dyngen)

G = c(1:50)
amt.genes = length(G)
tfs = G[1:20]
target2tfs <- lapply(G, function(g) {
  if(g%in%tfs) {
    return(double())
  } else {
    return(sample(tfs, 1+rbinom(1, length(tfs)-1, 0.02)))
  }
})
names(target2tfs) = G

# source("bin/ba_network.R")
# G = c(1:50)
# ba.network <- generate.ba(amnt.nodes = length(G), amnt.edges = 100)
# ba.net.df <- network.to.df(ba.network) # indien nodig
# tfs = G[G %in% ba.net.df$i]
# target2tfs = setNames(ba.network$neighbours, G)

## generating reaction formulas (the probability that a reaction occurs in [t, t dt])

kterms = c()
formulae <- unlist(recursive = F, lapply(G, function(g) {
  regs = target2tfs[[g]]#ba.network$neighbours[[g]]
  ## mRNA
  # production of mRNA
  if (length(regs) == 0) {
    mrnaprod.formula <- fvar("r", g) * fvar("a0", g)
  } else {
    inputs <- fsum(lapply(regs, function(r) fvar("y", r) / fvar("k", r, g)))
    numerator <- fvar("a0", g) + (fvar("a1") * inputs)
    denominator <- fcon(1) + inputs
    mrnaprod.formula <- fvar("r", g) * numerator / denominator
    
    kterms <<- c(kterms, lapply(regs, function(r) fvar("k", r, g)@string))
  }
  
  mrnaprod.nu <- rep(0, 2*length(G))
  mrnaprod.nu[[g]] <- 1
  
  # mRNA decay
  mrnadecay.formula <- fvar("d", g) * fvar("x", g)
  
  mrnadecay.nu <- rep(0, 2*length(G))
  mrnadecay.nu[[g]] <- -1
  
  ## Protein
  # production of protein
  protprod.formula <- fvar("p", g) * fvar("x", g)
  
  protprod.nu <- rep(0, 2*length(G))
  protprod.nu[[length(G)+g]] <- 1
  
  # production of protein
  protdecay.formula <- fvar("q", g) * fvar("y", g)
  
  protdecay.nu <- rep(0, 2*length(G))
  protdecay.nu[[length(G)+g]] <- -1
  
  # return formulae
  list(
    list(formula = mrnaprod.formula, nu = mrnaprod.nu),
    list(formula = mrnadecay.formula, nu = mrnadecay.nu),
    list(formula = protprod.formula, nu = protprod.nu),
    list(formula = protdecay.formula, nu = protdecay.nu)
  )
}))

formulae.strings <- sapply(formulae, function(fl) fl$formula@string)
formulae.nus <- sapply(formulae, function(fl) fl$nu)

#unlist(sapply(formulae, function(fl) {lapply(extract.variables(fl$formula), function(x) {x@string})}))

# 
# kterms = c()
# production = mapply(function(g, tfs) {
#   if (length(tfs) > 0) {
#     inputs = sapply(tfs, function(tf) paste0("(x", tf, "/k", tf, "g", g, ")"))
# 
#     kterms <<- c(kterms, sapply(tfs, function(tf) paste0("k", tf, "g", g)))
# 
#     up = paste0("(a0g", g, "+", paste("a1", inputs, sep="*", collapse="+"), ")")
#     down = paste0("(1+", paste0(inputs, collapse="+"), ")")
#   } else {
#     up = paste0("(a0g", g, ")")
#     down = "1"
#   }
#   return(paste0("r", g, "*", up, "/", down))
# }, tfs=target2tfs, g=names(target2tfs))
# 
# decay = mapply(function(g, tfs) {
#   return(paste0("d", g, "*x" , g))
# }, tfs=target2tfs, g=names(target2tfs))
#
#formulae = c(production, decay)

## generating the (initial) parameters of the system

R = rlnorm(length(G), log(10))/5
names(R) = paste0("r_", G)
D = rep(1,length(G))/5
names(D) = paste0("d_", G)
K = sapply(kterms, function(kterm) R[paste0("r_",str_replace(kterm, "k_(\\d*)_\\d*", "\\1"))]/2)
names(K) = kterms
A0 = rep(0.1, length(G))
names(A0) = paste0("a0_", G)

X0 = c(round(R/2, 0), rep(0, length(G)))
names(X0) = c(paste0("x_", G), paste0("y_", G))

P = rep(1,length(G))/10
names(P) = paste0("p_", G)
Q = rep(1,length(G))/10
names(Q) = paste0("q_", G)

## simulation

# function to postprocess the ssa output
process_ssa <- function(out, starttime=0) {
  final.time <- out$args$final.time
  data <- out$timeseries
  lastrow <- tail(data, n=1)
  lastrow$t <- final.time
  data[nrow(data)+1,] <- lastrow
  data <- data[data$t<=final.time,]
  times <- data$t + starttime
  data <- as.matrix(data[,-1])

  return(list(times=times, expression=data))
}

## Simulate multiple single cells following the same trajectory

# we will change the basal expression level of transcription factors during the simulation
# we need to determine these perturbations beforehand as they will be the same for every cell
# initial basal expression levels
totaltime = 40
burntime = 4

A0[] = 0.05
#A0[sample(tfs, 1+rbinom(1, length(tfs)-1, 0.2))] = 1
params = c(a1=1, R, D, K, A0, P, Q)
params = c(t=0, params)

# perturbations (saved in A0s)
# linear
params = matrix(params, nrow=nruns, byrow=T, ncol = length(params), dimnames = list(c(1:nruns), names(params)))
params[,"t"] = as.numeric(seq(0, totaltime-4, length=nrow(params)))
for(i in c(1:length(tfs))) {
  g = i
  start = i
  length = sample(seq(min(4, nrow(params)-start), nrow(params) - start), 1)
  for(i in seq(start, start+length)) {
    params[i,paste0("a0_", g)] = 1
  }
}

# add burn in time to params
params[,"t"] = params[,"t"] + burntime
params = rbind2(params[1,], params)
params[1,"t"] = 0

fit = function(X) apply(X, 2, function(x) (x-min(x))/(max(x)-min(x)))
pheatmap(fit(params[,apply(params, 2, sd) > 0]), cluster_cols=F, cluster_rows=F)

simulate_cell = function(timeofsampling=NULL) {
  requireNamespace("fastgssa")
  
  if (!is.null(timeofsampling)) {
    time = timeofsampling+burntime
  } else {
    time = totaltime+burntime
  }

  set.seed(1)
  out <- fastgssa::ssa(X0,formulae.strings,formulae.nus, burntime+time,params, method=fastgssa::ssa.direct(), recalculate.all =F)
  output = process_ssa(out)
  expression = output$expression
  times = output$times
  burnin = times > burntime
  expression = output$expression[burnin,]
  times = times[burnin] - burntime

  # return either the full expression matrix, or the expression at timeofsampling
  if(is.null(timeofsampling)) {
    rownames(expression) = c(1:nrow(expression))
    return(list(expression=expression, times=times))
  } else {
    #if(timeofsampling > last(times)) warning("samplingtime larger than simulated time")
    return(tail(expression, n=1)[1,])
    #return(expression[max(1, findInterval(timeofsampling, times)),])
  }
}

celltimes = runif(200, 0, totaltime)

#library(profvis)
#profvis({simulate_cell(1)})

simulate_cell(40)

cells = mclapply(celltimes, simulate_cell, mc.cores=1)
cells = qsub.lapply(celltimes, simulate_cell)

expression = matrix(unlist(cells), nrow=length(cells), byrow=T, dimnames = list(c(1:length(cells)), names(cells[[1]])))
E = expression[,str_detect(colnames(expression), "x_")]
E = E[,apply(E, 2, sd) > 0]

Eprot = expression[,str_detect(colnames(expression), "y_")]
Eprot = Eprot[,apply(Eprot, 2, sd) > 0]

pheatmap(t(E), scale="row", cluster_cols=F, cluster_rows=T, clustering_distance_rows = "correlation", annotation_col = data.frame(time=celltimes, row.names = rownames(E)))
pheatmap(t(Eprot), scale="row", cluster_cols=F, cluster_rows=T, clustering_distance_rows = "correlation", annotation_col = data.frame(time=celltimes, row.names = rownames(E)))

library(ggplot2)
datadf = melt(E)
colnames(datadf) = c("time", "gene", "count")

datadf = datadf[datadf$gene %in% paste0("x_", tfs),]
ggplot(datadf) + geom_line(aes(time, count, group=gene, color=gene))

library(SCORPIUS)
space = reduce.dimensionality(correlation.distance(E),ndim = 3)
#
# space = tsne::tsne(as.dist(correlation.distance(E)))
# colnames(space) = c("Comp1", "Comp2")
#
# space = diffusionMap::diffuse(as.dist(correlation.distance(E)))
# space = space$X[,c(1:2)]
# colnames(space) = c("Comp1", "Comp2")

trajectory = infer.trajectory(space)
draw.trajectory.plot(space, celltimes, trajectory$final.path) + scale_colour_distiller(palette = "RdYlBu")
rownames(E) = c(1:nrow(E))
draw.trajectory.heatmap(E, trajectory$time, as.factor(cut(celltimes, breaks=length(celltimes), labels=F)))

SCORPIUS::evaluate.trajectory(trajectory$time, celltimes)

draw.trajectory.heatmap(E, celltimes, factor(cut(celltimes, breaks=length(celltimes), labels=F)), show.labels.row = T)


# check correlations between cells
cellcor = cor(t(E[order(celltimes),]))
pheatmap(cellcor, cluster_rows=F, cluster_cols=F)

# check correlation between perturbations and simulated expression
runtimes = seq(0, totaltime, length=nruns+1)
Eperturb = matrix(unlist(sapply(celltimes, function(t) A0s[min(length(A0s),max(0, findInterval(t, runtimes)))])), nrow=length(celltimes), byrow = T)
colnames(Eperturb) = paste0("x", G)

E = matrix(unlist(cells), nrow=length(cells), byrow=T, dimnames = list(c(1:length(cells)), names(cells[[1]])))

joined = left_join(melt(E, value.name="simul"), melt(Eperturb, value.name = "perturb"), c("Var1", "Var2"))
colnames(joined)[1:2] = c("time", "gene")

ggplot(joined) + geom_point(aes(perturb, simul)) + facet_wrap(~gene, scales = c("free_y"))

pheatmap(fit(params[,apply(params, 2, sd) > 0]), cluster_cols=F, cluster_rows=F)

# check a random tf and his targets, for example to check the time lag
tf = tfs[[1]]
targets = G[sapply(target2tfs, function(tfs) (tf %in% tfs) && length(tfs)==1)]
length(targets)
pheatmap(t(E[order(celltimes),c(tf,targets)]), scale="row", cluster_rows=F, cluster_cols = F)

# check protein and mRNA
i = 23
x = paste0("x_",i)
y = paste0("y_",i)
plotdata = bind_rows(
  data.frame(expression=E[order(celltimes), x], g=x, time=sort(celltimes)),
  data.frame(expression=Eprot[order(celltimes), y], g=y, time=sort(celltimes))
)

ggplot(plotdata) + geom_line(aes(time, expression, group=g, color=g))
