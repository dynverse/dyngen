.ba.initialise.network <- function(amnt.nodes) {
  degree <- rep(0, amnt.nodes)
  forw.neighbours <- lapply(seq_len(amnt.nodes), function(i) numeric(0))
  back.neighbours <- lapply(seq_len(amnt.nodes), function(i) numeric(0))
  list(degree = degree, neighbours = forw.neighbours, back.neighbours = back.neighbours)
}

.ba.add.edge <- function(net, i, j) {
  net$degree[c(i,j)] <- net$degree[c(i,j)] + 1
  net$neighbours[[i]] <- c(net$neighbours[[i]], j)
  net$back.neighbours[[j]] <- c(net$back.neighbours[[j]], i)
  net
}

.ba.remove.edge <- function(net, i, j) {
  net$degree[c(i,j)] <- net$degree[c(i,j)] - 1
  n <- net$neighbours[[i]]
  n <- n[n != j]
  net$neighbours[[i]] <- n
  m <- net$back.neighbours[[j]]
  m <- m[m != i]
  net$back.neighbours[[j]] <- m
  net
}

.ba.remove.node <- function(net, i) {
  n <- net$neighbours[[i]]
  for (j in n) {
    net <- .ba.remove.edge(net, i, j)
  }
  net
}

.ba.add.node.with.edges <- function(net, i, num.edges, offset.exponent) {
  if (num.edges >= i) stop(sQuote("num.edges"), " must be smaller than ", sQuote("i"))
  
  bn <- net$back.neighbours[[i]]
  poss.neighs <- seq_len(i-1) # only smaller than i
  if (length(bn) > 0) {
    poss.neighs <- poss.neighs[-bn]
  }
  w <- net$degree[poss.neighs]
  w[w == 0] <- 1
  w <- (w / sum(w)) 
  w <- w ^ offset.exponent
  neighs <- sample(poss.neighs, size = num.edges, replace=F, prob=w)
  for (j in neighs) {
    net <- .ba.add.edge(net, i, j)
  }
  net
}

# adapted from http://stackoverflow.com/questions/24845909/generate-n-random-integers-that-sum-to-m-in-r
.ba.num.edges <- function(amnt.nodes, amnt.edges, sd = 1, pos.only = TRUE) {
  vec <- rnorm(amnt.nodes, amnt.edges/amnt.nodes, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * amnt.edges)
  deviation <- amnt.edges - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(amnt.nodes, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  for (i in seq_along(vec)) {
    if (vec[[i]] >= i) {
      num.edges.to.redistribute <- vec[[i]] - i + 1
      vec[[i]] <- i - 1
      places.to.put <- sample(seq(i+1, length(vec)), size = num.edges.to.redistribute, replace = T)
      for (j in places.to.put) {
        vec[[j]] <- vec[[j]] + 1
      }
    }
  }
  vec
}

generate.ba <- function(amnt.nodes, amnt.edges, offset.exponent = 1, trace=T) {
  if (amnt.edges > amnt.nodes * (amnt.nodes - 1) / 2) stop(sQuote("amnt.edges"), " is too large, as a graph with N nodes can only contain N*(N-1)/2 edges")
  
  net <- .ba.initialise.network(amnt.nodes)
  
  # generate number of outgoing edges per node
  nums <- .ba.num.edges(amnt.nodes, amnt.edges)
  
  # add nodes
  for (i in seq_len(amnt.nodes)) {
    if (nums[[i]] > 0) {
      net <- .ba.add.node.with.edges(net, i, nums[[i]], offset.exponent)
    }
  }
  
  # transform network into a data frame
  network <- bind_rows(lapply(seq_len(amnt.nodes), function(i) {
    data.frame(i = net$neighbours[[i]], j = rep(i, length(net$neighbours[[i]])))
  }))
  
  network
}
