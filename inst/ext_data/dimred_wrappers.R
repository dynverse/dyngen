dimred_pca = function(x, ndim=3) {
  space = prcomp(t(x))$rotation[,seq_len(ndim)]
  process_dimred(space)
}
dimred_simlr = function(x, ndim=3, nclusters=4) {
  result = SIMLR::SIMLR(t(x), nclusters)
  S = result$S
  space = tsne::tsne(as.dist(max(S)-S), k = ndim)
  process_dimred(space)
}

dimred_mds = function(x, ndim=3) {
  space = SCORPIUS::reduce_dimensionality(SCORPIUS::correlation_distance(x),ndim = ndim)
  process_dimred(space)
}

dimred_mds_sammon = function(x, ndim=3) {
  dist = SCORPIUS::correlation_distance(x)
  space <- MASS::sammon(dist, k = ndim)$points
  process_dimred(space)
}

dimred_mds_isomds = function(x, ndim=3) {
  dist = SCORPIUS::correlation_distance(x)
  space <- MASS::isoMDS(dist, k = ndim)$points
  process_dimred(space)
}

#dimred_lmds = function(x, ndim=3) {
#  mds.out <- dambiutils::mds_withlandmarks(x %>% as.data.frame, SCORPIUS::correlation_distance, k = ndim, landmark.method = "naive", num.landmarks = min(1000, round(nrow(x)*0.1)), num.seed.landmarks = 10, pca.normalisation = FALSE)
#  process_dimred(mds.out$S)
#}

dimred_mds_smacof = function(x, ndim=3) {
  dist = SCORPIUS::correlation_distance(x)
  space <- smacof::mds(as.dist(dist), type = "ratio", ndim = ndim)$conf
  process_dimred(space)
}

dimred_tsne = function(x, ndim=3) {
  space = Rtsne::Rtsne(as.dist(SCORPIUS::correlation_distance(x)), dims = ndim, is_distance = TRUE,  perplexity = 20)$Y
  rownames(space) = rownames(x)
  process_dimred(space)
}

dimred_dp = function(x, ndim=3, neigen=ndim) {
  space = diffusionMap::diffuse(as.dist(SCORPIUS::correlation_distance(x)), neigen=neigen, maxdim=ndim, delta=10e-3)
  process_dimred(space$X[,seq_len(ndim)])
}

dimred_ica = function(x, ndim=3) {
  space = fastICA::fastICA(t(scale(t(x))), ndim)$S
  process_dimred(space)
}

dimred_lle = function(x, ndim=3) {
  k = lle::calc_k(t(scale(t(x))), ndim)
  k = k$k[which.min(k$rho)]
  space = lle::lle(t(scale(t(x))), ndim, k)$Y
  process_dimred(space)
}

process_dimred = function(space) {
  space = as.matrix(space)
  colnames(space) = paste0("Comp", 1:ncol(space))
  space
}


