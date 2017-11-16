## TODO:
# * doublets, using an extracellrate
#     What is the cell state & other gold standard values?


## scRNAseq functions
#' @importFrom stats rbinom
amplify = function(a=0, steps=100, rate=0.05) {
  for (i in 1:steps) {
    a = a + stats::rbinom(length(a), a, rate)
  }
  a
}

#' @importFrom stats rmultinom runif
libprep = function(counts, lysisrate = 0.6, capturerate = 0.1, amplify = T, amplifysteps = 100, amplifyrate = c(0.0001, 0.03), sequencerate=0.1, cellcapturerate = 1, celldoubletrate = 0, verbose=F, verbose_plot_cell=50, verbose_follow_gene="x_1") {
  curcounts = counts
  curcounts = counts_cellcaptured = counts[stats::runif(nrow(counts)) <= cellcapturerate, ]
  
  if (verbose && !(verbose_plot_cell %in% rownames(curcounts))) {verbose_plot_cell = rownames(curcounts)[[1]]}
  
  print(curcounts %>% dim)
  
  curcounts = counts_lysed = apply(curcounts, 1, function(col) stats::rmultinom(1, sum(col)*lysisrate, col)[,1]) %>% t
  dimnames(curcounts) = dimnames(counts_lysed) = dimnames(counts_cellcaptured)
  curcounts = counts_captured = apply(curcounts, 1, function(col) stats::rmultinom(1, sum(col)*capturerate, col)[,1]) %>% t
  dimnames(curcounts) = dimnames(counts_captured) = dimnames(counts_cellcaptured)
  
  # different rates per cell and per gene
  amplifyrates_cells = stats::runif(nrow(curcounts), min = amplifyrate[1], max=amplifyrate[2])
  amplifyrates_genes = stats::runif(ncol(curcounts), min = amplifyrate[1], max=amplifyrate[2])
  amplifyrates = sqrt(amplifyrates_cells %o% amplifyrates_genes)
  
  if(amplify) {
    curcounts = counts_amplified = lapply(1:nrow(curcounts), function(row_id) {
      amplify(curcounts[row_id, ], amplifysteps, amplifyrates[row_id, ])
    }) %>% do.call(rbind, .)
    dimnames(curcounts) = dimnames(counts_amplified) = dimnames(counts_cellcaptured)
  }
  
  curcounts = counts_sequences = apply(curcounts, 1, function(col) stats::rmultinom(1, sum(col)*sequencerate, col+0.000001)[,1]) %>% t # pseudocounts added because sum of probs cannot be zero
  dimnames(curcounts) = dimnames(counts_sequences) = dimnames(counts_cellcaptured)
  
  if(verbose) {
    overview = tibble()
    overview = bind_rows(
      tibble(gene=colnames(counts_cellcaptured), count = counts_cellcaptured[verbose_plot_cell, ], step="cell"),
      tibble(gene=colnames(counts_cellcaptured), count = counts_lysed[verbose_plot_cell, ], step="lysed"),
      tibble(gene=colnames(counts_cellcaptured), count = counts_captured[verbose_plot_cell, ], step="captured"),
      tibble(gene=colnames(counts_cellcaptured), count = counts_amplified[verbose_plot_cell, ], step="amplified"),
      tibble(gene=colnames(counts_cellcaptured), count = counts_sequences[verbose_plot_cell, ], step="sequenced")
    ) %>% mutate(step=factor(step, levels=unique(step)))
    plot = overview %>% 
      ggplot() +
      geom_histogram(aes(count), bins=50) + 
      facet_wrap(~step, ncol=1) + 
      geom_vline(aes(xintercept=count, color=gene), data=overview %>% filter(gene %in% verbose_follow_gene))
    print(plot)
  }
  curcounts
}

simulate_scrnaseq = function(counts, platform) {
  libprep(counts, amplifyrate = c(platform$amplifyrate_min, platform$amplifyrate_max), lysisrate=platform$lysisrate, capturerate=platform$capturerate, amplify=platform$amplify, cellcapturerate=platform$cellcapturerate, verbose=F, verbose_plot_cell= "C1", verbose_follow_gene=c("G3","G10", "G50"))
}
