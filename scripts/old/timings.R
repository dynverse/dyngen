library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(pheatmap)

G = c(1:50)
tfs = G[1:10]
target2tfs <- lapply(G, function(g) sample(tfs, 1+rbinom(1, length(tfs)-1, 0.1)))
names(target2tfs) = G

interactionid = 0
kterms = c()
production = mapply(function(g, tfs) {
  inputs = sapply(tfs, function(tf) paste0("(x", tf, "/k", tf, "g", g, ")"))

  kterms <<- c(kterms, sapply(tfs, function(tf) paste0("k", tf, "g", g)))

  up = paste0("(a0g", g, "+", paste("a1", inputs, sep="*", collapse="+"), ")")
  down = paste0("(1+", paste0(inputs, collapse="+"), ")")

  #up = paste0(sapply(tfs, function(tf) paste0("(X", tf, "/k", tf, "g", g, ")")), collapse="+")
  return(paste0("r", g, "*", up, "/", down))
}, tfs=target2tfs, g=names(target2tfs))

decay = mapply(function(g, tfs) {
  return(paste0("d", g, "*x" , g))
}, tfs=target2tfs, g=names(target2tfs))

formulas = c(production, decay)

R = rlnorm(length(G), log(10))
names(R) = paste0("r", G)
D = rep(1,length(G))
names(D) = paste0("d", G)
K = sapply(kterms, function(kterm) R[paste0("r",str_replace("k1g2", "k(\\d*)g\\d*", "\\1"))])
names(K) = kterms
A0 = rep(0.1, length(G))
names(A0) = paste0("a0g", G)
A0[sample(tfs, 1+rbinom(1, length(tfs)-1, 0.4))] = 1

params = c(a0=0.1, a1=1, R, D, K, A0)
nu <- cbind2(diag(nrow=length(G), ncol=length(G)), -diag(nrow=length(G), ncol=length(G)))
X0 = round(R/4, 0)
names(X0) = paste0("x", G)

tf = 10

# x0 <- X0
# a <- formulas
# parms <- params
# method <- "OTL"
# verbose <- T
# consoleInterval <- 100

out <- fastgssa::ssa(x0 = X0, a = formulas, nu = nu, parms = params, tf = tf, method="OTL", verbose = T, consoleInterval = 100)
out$stats$elapsedWallTime %>% as.numeric()

out <- GillespieSSA::ssa(x0 = X0, a = formulas, nu = nu, parms = params, tf = tf, method="OTL", verbose = T, consoleInterval = 100)
out$stats$elapsedWallTime[["elapsed"]]

library(pbapply)

fastgssas <- pbsapply(seq_len(10), function(i) {
  out <- fastgssa::ssa(x0 = X0, a = formulas, nu = nu, parms = params, tf = tf, method="OTL", verbose = F)
  out$stats$elapsedWallTime %>% as.numeric()
})
gillespies <- pbsapply(seq_len(10), function(i) {
  out <- GillespieSSA::ssa(x0 = X0, a = formulas, nu = nu, parms = params, tf = tf, method="OTL", verbose = F)
  out$stats$elapsedWallTime[["elapsed"]]
})
timings <- bind_rows(data.frame(package = "fastgssa", time = fastgssas), data.frame(package = "GillespieSSA", time = gillespies))
ggplot(timings, aes(package, time)) + geom_boxplot() + theme_classic()
