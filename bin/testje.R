library(GillespieSSA)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(pheatmap)

G = c(1:50)
tfs = G[1:10]
target2tfs <- lapply(G, function(g) {
  if(g%in%tfs) {
    return(double())
  } else {
    return(sample(tfs, 1+rbinom(1, length(tfs)-1, 0.1)))
  }
})
names(target2tfs) = G

# source("ba_network.R")
# ba.network <- generate.ba(amnt.nodes = 10, amnt.edges = 20)
# ba.net.df <- network.to.df(ba.network) # indien nodig

## generating reaction formulas (the probability that a reaction occurs in [t, t dt])

interactionid = 0
kterms = c()
production = mapply(function(g, tfs) {
  if (length(tfs) > 0) {
    inputs = sapply(tfs, function(tf) paste0("(x", tf, "/k", tf, "g", g, ")"))
    
    kterms <<- c(kterms, sapply(tfs, function(tf) paste0("k", tf, "g", g)))
    
    up = paste0("(a0g", g, "+", paste("a1", inputs, sep="*", collapse="+"), ")")
    down = paste0("(1+", paste0(inputs, collapse="+"), ")")
  } else {
    up = paste0("(a0g", g, ")")
    down = "1"
  }
  
  #up = paste0(sapply(tfs, function(tf) paste0("(X", tf, "/k", tf, "g", g, ")")), collapse="+")
  return(paste0("r", g, "*", up, "/", down))
}, tfs=target2tfs, g=names(target2tfs))

decay = mapply(function(g, tfs) {
  return(paste0("d", g, "*x" , g))
}, tfs=target2tfs, g=names(target2tfs))

formulas = c(production, decay)

## generating the (initial) parameters of the system

R = rlnorm(length(G), log(10))
names(R) = paste0("r", G)
D = rep(1,length(G))
names(D) = paste0("d", G)
K = sapply(kterms, function(kterm) R[paste0("r",str_replace(kterm, "k(\\d*)g\\d*", "\\1"))]/2)
names(K) = kterms
A0 = rep(0.1, length(G))
names(A0) = paste0("a0g", G)

nu <- cbind2(diag(nrow=length(G), ncol=length(G)), -diag(nrow=length(G), ncol=length(G)))
X0 = round(R/2, 0)
names(X0) = paste0("x", G)

## simulation

tf = 5

# function to postprocess the ssa output
process_ssa = function(out, starttime=0) {
  tf = out$args$tf
  rownames(out$data) = NULL
  data = as.data.frame(out$data)
  lastrow = tail(data, n=1)
  lastrow$V1 = tf
  data[nrow(data)+1,] <- lastrow
  data <- data[data$V1<tf,]
  times = data$V1 + starttime
  data = as.matrix(data[,-1])
  
  rownames(data) = NULL
  
  return(list(times=times, expression=data))
}

# we need to perturb the system during differentiation, but ssa does not allow this
# we therefore combine different ssa runs and each time perturb parameters and set X0 to Xn of the previous run

# first run
A0[] = 0.05
A0[sample(tfs, 1+rbinom(1, length(tfs)-1, 0.5))] = 1
params = c(a1=1, R, D, K, A0)

print(A0)

runtimes = c()

tf = 5
burntf = 4
out <- ssa(X0,formulas,nu,params,tf=tf, method = "BTL")
output = process_ssa(out)
expression = output$expression
times = output$times
burnin = times > burntf
expression = output$expression[burnin,]
times = times[burnin]
runtimes[length(runtimes)+1] = last(times)
runs = rep(1, length(times))

# other runs: perturb TF activity
tf=1
lapply(c(1:20), function(i) {
  X0 = tail(expression, n=1)[1,]
  
  A0[sample(tfs, 1)] <<- 1
  A0[sample(tfs, 1)] <<- 0.05
  #A0[sample(tfs, 1+rbinom(1, length(tfs)-1, 0.5))] = 1
  print(A0)
  params = c(a1=1, R, D, K, A0)
  out <- ssa(X0,formulas,nu,params,tf=tf, method="BTL")
  output = process_ssa(out, last(times))
  expression<<-rbind2(expression, output$expression)
  times <<- c(times, output$times)
  runtimes[length(runtimes)+1] <<- last(times)
  runs <<- c(runs, rep(i+1, length(output$times)))
  return("yay!")
})

rownames(expression) = c(1:nrow(expression))
pheatmap(t(expression), scale="row", cluster_cols=F, cluster_rows=T, clustering_distance_rows = "correlation", annotation_col = data.frame(run=factor(runs), row.names = rownames(expression)))

library(ggplot2)
datadf = melt(expression)
colnames(datadf) = c("sample", "gene", "count")

datadf = datadf[datadf$gene %in% c("x1", "x2", "x3", "x4", "x5"),]

ggplot(datadf) + geom_line(aes(sample, count, group=gene, color=gene))


library(SCORPIUS)

samples = sample(c(1:nrow(expression)), 500, replace = F)
sampletimes = times[samples]

E = expression[samples, ]

space = reduce.dimensionality(correlation.distance(E),ndim = 2)
trajectory = infer.trajectory(space)
draw.trajectory.plot(space, sampletimes, trajectory$final.path) + scale_colour_distiller(palette = "RdYlBu")
rownames(E) = c(1:nrow(E))
draw.trajectory.heatmap(E, trajectory$time, as.factor(cut(sampletimes, breaks=99, labels=F)))

draw.trajectory.heatmap(E, sampletimes, as.factor(cut(sampletimes, breaks=99, labels=F)), show.labels.row = T)
