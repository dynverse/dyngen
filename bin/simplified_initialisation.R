library(fastgssa)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(pheatmap)
library(parallel)
library(PRISM)

amnt.genes <- 50
amnt.edges <- 200
ba.network <- generate.ba(amnt.nodes = amnt.genes, amnt.edges = amnt.edges, reverse.edges = T, offset.exponent = 1.5)


formulae <- sapply(seq_len(amnt.genes), function(g) {
  list(
    generate.production.formula(g, ba.network$incoming.edges[[g]], amnt.genes),
    generate.decay.formula(g, amnt.genes)
  )
})

formulae.vars <- lapply(formulae, function(fl) fl$variables)
formulae.strings <- sapply(formulae, function(fl) fl$formula)
formulae.nus <- sapply(formulae, function(fl) fl$nu)



## generating the (initial) parameters of the system
gene.ix <- seq_len(amnt.genes)
R = setNames(rlnorm(length(G), log(10)), paste0("r", gene.ix))
D = setNames(rep(1,length(G)), paste0("d", gene.ix))

K = sapply(kterms, function(kterm) R[paste0("r",str_replace(kterm, "k(\\d*)g\\d*", "\\1"))]/2)
names(K) = kterms
A0 = rep(0.1, length(G))
names(A0) = paste0("a0g", G)

nu <- cbind2(diag(nrow=length(G), ncol=length(G)), -diag(nrow=length(G), ncol=length(G)))
X0 = round(R/2, 0)
names(X0) = paste0("x", G)

