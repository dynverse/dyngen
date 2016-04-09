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

G = c(1,2,3)
tfs = c(1,2)
target2tfs = list(numeric(), c(1), c(2))

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

R = rlnorm(length(G), log(20))/5
R[] = 20
names(R) = paste0("r_", G)
D = rep(1,length(G))/5
names(D) = paste0("d_", G)
K = sapply(kterms, function(kterm) (R/D)[paste0("r_",str_replace(kterm, "k_\\d*_(\\d*)", "\\1"))]/4)
names(K) = kterms
A0 = rep(0.05, length(G))
names(A0) = paste0("a0_", G)
A0[[1]] = 1

X0 = c(rep(0, length(G)), rep(0, length(G)))
X0[[1]] = 1
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
burntime = 0

params = c(a1=1, R, D, K, A0, P, Q)
params = c(t=0, params)
params = matrix(params, nrow=1, byrow=T, ncol = length(params), dimnames = list(1, names(params)))

simulate_cell = function(timeofsampling=NULL) {
  requireNamespace("fastgssa")
  
  if (!is.null(timeofsampling)) {
    time = timeofsampling
  } else {
    time = totaltime
  }
  
  out <- fastgssa::ssa(X0,formulae.strings,formulae.nus, burntime+time, params, method=fastgssa::ssa.direct(), recalculate.all =F)
  output = process_ssa(out)
  expression = output$expression
  times = output$times
  burnin = times > burntime
  expression = methods::rbind2(output$expression[max(1, findInterval(burntime, output$times)),], output$expression[burnin,])
  times = c(0,times[burnin] - burntime)
  
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

celltimes = runif(100, 0, totaltime)

cells = mclapply(celltimes, simulate_cell, mc.cores=1)
cells = qsub.lapply(celltimes, simulate_cell)

expression = matrix(unlist(cells), nrow=length(cells), byrow=T, dimnames = list(c(1:length(cells)), names(cells[[1]])))
E = expression[,str_detect(colnames(expression), "x_")]
E = E[,apply(E, 2, sd) > 0]

Eprot = expression[,str_detect(colnames(expression), "y_")]
Eprot = Eprot[,apply(Eprot, 2, sd) > 0]

pheatmap(t(E[order(celltimes),]), scale="row", cluster_cols=F, cluster_rows=T, clustering_distance_rows = "correlation", annotation_col = data.frame(time=celltimes, row.names = rownames(E)))
pheatmap(t(Eprot[order(celltimes),]), scale="row", cluster_cols=F, cluster_rows=T, clustering_distance_rows = "correlation", annotation_col = data.frame(time=celltimes, row.names = rownames(E)))

library(ggplot2)
datadf = melt(Eprot[order(celltimes),])
colnames(datadf) = c("time", "gene", "count")

ggplot(datadf) + geom_line(aes(time, count, group=gene, color=gene))

# check a random tf and his targets, for example to check the time lag
tf = tfs[[1]]
targets = G[sapply(target2tfs, function(tfs) (tf %in% tfs) && length(tfs)==1)]
length(targets)
pheatmap(t(E[order(celltimes),c(tf,targets)]), scale="row", cluster_rows=F, cluster_cols = F)

# check protein and mRNA
i = 3
x = paste0("x_",i)
y = paste0("y_",i)
plotdata = bind_rows(
  data.frame(expression=E[order(celltimes), x], g=x, time=sort(celltimes)),
  data.frame(expression=Eprot[order(celltimes), y], g=y, time=sort(celltimes))
)

ggplot(plotdata) + geom_point(aes(time, expression, group=g, color=g))

#

nknots = 8
Esmooth = cbind2(apply(E, 2, function(x) {smooth.spline(celltimes, x, nknots=nknots)$y}), apply(Eprot, 2, function(x) {smooth.spline(celltimes,x, nknots=nknots)$y}))
Esmooth = cbind2(apply(E[order(celltimes), ], 2, smooth, kind="3RS3R"), apply(Eprot[order(celltimes), ], 2, smooth, kind="3RS3R"))
plotdata = melt(Esmooth, value.name="expression")
plotdata = melt(cbind2(E, Eprot), value.name="expression")
colnames(plotdata) = c("time", "var", "expression")
plotdata$time = celltimes[plotdata$time]
plotdata$type = sapply(str_split(plotdata$var, "_"), function(x) x[[1]])
plotdata$gene = sapply(str_split(plotdata$var, "_"), function(x) x[[2]])

ggplot(plotdata, aes(x=time, y=expression)) + 
  geom_point(aes(color=gene, shape=type), alpha=1)+
  geom_smooth(aes(group=var, color=gene, linetype=type)) + 
  facet_wrap(~gene)

# check multiple individual cells
plotdata = bind_rows(lapply(seq_len(10), function(i) {
  out = simulate_cell()
  data.frame(expression=out$expression[,x], simulation=i, time=out$times)
}))
ggplot(plotdata) + geom_line(aes(time, expression, group=simulation, color=simulation))
