dimreds = list("pca"=pca, "mds"=mds, "tsne"=tsne, "dp"=dp, "ica"=ica, "lle"=lle)

pca = function(x, ndim=3) {
  space = prcomp(t(x))$rotation[,seq_len(ndim)]
  colnames(space) = paste0("Comp", 1:ncol(space))
  space
}

mds = function(x, ndim=3) {
  SCORPIUS::reduce.dimensionality(SCORPIUS::correlation.distance(x),ndim = ndim)
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
  fastICA::fastICA(t(scale(t(x))), ndim)$S
}

lle = function(x, ndim=3) {
  k = lle::calc_k(t(scale(t(x))), ndim)
  k = k$k[which.min(k$rho)]
  space = lle::lle(t(scale(t(x))), ndim, k)$Y
  space
}