## scRNAseq functions
amplify = function(a=0, steps=100, rate=0.05) {
  for (i in 1:steps) {
    a = a + rbinom(length(a), a, rate)
  }
  a
}

libprep = function(counts, lysisrate = 0.6, capturerate = 0.1, amplify = T, amplifysteps = 100, amplifyrate = c(0.0001, 0.03), sequencerate=0.1, verbose=F, verbose_plot_cell=50, verbose_follow_gene="x_1") {
  curcounts = counts
  curcounts = counts_lysed = apply(curcounts, 1, function(col) rmultinom(1, sum(col)*lysisrate, col)[,1]) %>% t
  dimnames(curcounts) = dimnames(counts_lysed) = dimnames(counts)
  curcounts = counts_captured = apply(curcounts, 1, function(col) rmultinom(1, sum(col)*capturerate, col)[,1]) %>% t
  dimnames(curcounts) = dimnames(counts_captured) = dimnames(counts)
  
  # different rates per cell and per gene
  amplifyrates_cells = runif(nrow(curcounts), min = amplifyrate[1], max=amplifyrate[2])
  amplifyrates_genes = runif(ncol(curcounts), min = amplifyrate[1], max=amplifyrate[2])
  amplifyrates = sqrt(amplifyrates_cells %o% amplifyrates_genes)
  
  if(amplify) {
    curcounts = counts_amplified = lapply(1:nrow(curcounts), function(rowid) {
      amplify(curcounts[rowid, ], amplifysteps, amplifyrates[rowid, ])
    }) %>% do.call(rbind, .)
  }
  dimnames(curcounts) = dimnames(counts_amplified) = dimnames(counts)
  
  curcounts = counts_sequences = apply(curcounts, 1, function(col) rmultinom(1, sum(col)*sequencerate, col+0.000001)[,1]) %>% t # pseudocounts added because sum of probs cannot be zero
  dimnames(curcounts) = dimnames(counts_sequences) = dimnames(counts)
  
  if(verbose) {
    overview = tibble()
    overview = bind_rows(
      tibble(gene=colnames(counts), count = counts[verbose_plot_cell, ], step="cell"),
      tibble(gene=colnames(counts), count = counts_lysed[verbose_plot_cell, ], step="lysed"),
      tibble(gene=colnames(counts), count = counts_captured[verbose_plot_cell, ], step="captured"),
      tibble(gene=colnames(counts), count = counts_amplified[verbose_plot_cell, ], step="amplified"),
      tibble(gene=colnames(counts), count = counts_sequences[verbose_plot_cell, ], step="sequenced")
    ) %>% mutate(step=factor(step, levels=unique(step)))
    plot = overview %>% ggplot() + geom_histogram(aes(count), bins=50) + facet_wrap(~step, ncol=1) + geom_vline(aes(xintercept=count, color=gene), data=overview %>% filter(gene %in% verbose_follow_gene))
    plot(plot)
  }
  curcounts
}

simulate_scrnaseq = function(expression, platform) {
  cellcounts = expression %>% {.*100} %>% round() %>% abs()
  libprep(cellcounts, amplifyrate = c(platform$amplifyrate_min, platform$amplifyrate_max), lysisrate=platform$lysisrate, capturerate=platform$capturerate, amplify=platform$amplify, verbose=T, verbose_plot_cell= "C1", verbose_follow_gene=c("G3","G10", "G50"))
}
