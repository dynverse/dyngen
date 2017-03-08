dimreds = list("pca"=pca, "mds"=mds, "tsne"=tsne, "dp"=dp, "ica"=ica, "lle"=lle)
dimreds = list("mds"=mds, "pca"=pca, "ica"=ica)

pca = function(x, ndim=3) {
  space = prcomp(t(x))$rotation[,seq_len(ndim)]
  colnames(space) = paste0("Comp", 1:ncol(space))
  space
}

mds = function(x, ndim=3) {
  space = SCORPIUS::reduce.dimensionality(SCORPIUS::correlation.distance(x),ndim = ndim)
  colnames(space) = paste0("Comp", 1:ncol(space))
  space
}

mds_sammon = function(x, ndim=3) {
  dist = SCORPIUS::correlation.distance(x)
  space <- MASS::sammon(dist, k = ndim)$points
  colnames(space) = paste0("Comp", 1:ncol(space))
  space
}

mds_isomds = function(x, ndim=3) {
  dist = SCORPIUS::correlation.distance(x)
  space <- MASS::isoMDS(dist, k = ndim)$points
  colnames(space) = paste0("Comp", 1:ncol(space))
  space
}

lmds = function(x, ndim=3) {
  mds.out <- dambiutils::mds_withlandmarks(x %>% as.data.frame, SCORPIUS::correlation.distance, k = ndim, landmark.method = "naive", num.landmarks = min(1000, round(nrow(x)*0.1)), num.seed.landmarks = 10, pca.normalisation = F)
  space <- mds.out$S
  colnames(space) = paste0("Comp", 1:ncol(space))
  space
}

mds_smacof = function(x, ndim=3) {
  dist = SCORPIUS::correlation.distance(x)
  space <- smacof::mds(as.dist(dist), type = "ratio", ndim = ndim)$conf
  colnames(space) = paste0("Comp", 1:ncol(space))
  space
}

tsne = function(x, ndim=3) {
  tsne::tsne(as.dist(SCORPIUS::correlation.distance(x)), k = ndim)
}

dp = function(x, ndim=3, neigen=3) {
  space = diffusionMap::diffuse(as.dist(SCORPIUS::correlation.distance(x)), neigen=neigen)
  space = space$X[,seq_len(ndim)]
  space
}

ica = function(x, ndim=3) {
  space = fastICA::fastICA(t(scale(t(x))), ndim)$S
  
  colnames(space) = paste0("Comp", 1:ncol(space))
  space
}

lle = function(x, ndim=3) {
  k = lle::calc_k(t(scale(t(x))), ndim)
  k = k$k[which.min(k$rho)]
  space = lle::lle(t(scale(t(x))), ndim, k)$Y
  space
}