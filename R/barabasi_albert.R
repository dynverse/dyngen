.ba.initialise.network <- function(amnt.nodes) {
  degree <- rep(0, amnt.nodes)
  in.degree <- rep(0, amnt.nodes)
  out.degree <- rep(0, amnt.nodes)
  incoming.edges <- lapply(seq_len(amnt.nodes), function(i) numeric(0))
  outgoing.edges <- lapply(seq_len(amnt.nodes), function(i) numeric(0))
  list(
    amnt.nodes = amnt.nodes, 
    degree = degree, 
    in.degree = in.degree, 
    out.degree = out.degree, 
    incoming.edges = incoming.edges, 
    outgoing.edges = outgoing.edges
  )
}

.ba.add.edge <- function(net, i, j) {
  net$degree[c(i,j)] <- net$degree[c(i,j)] + 1
  net$out.degree[[i]] <- net$out.degree[[i]] + 1
  net$outgoing.edges[[i]] <- c(net$outgoing.edges[[i]], j)
  net$in.degree[[j]] <- net$in.degree[[j]] + 1
  net$incoming.edges[[j]] <- c(net$incoming.edges[[j]], i)
  net
}

.ba.remove.edge <- function(net, i, j) {
  net$degree[c(i,j)] <- net$degree[c(i,j)] - 1
  net$out.degree[[i]] <- net$out.degree[[i]] - 1
  net$in.degree[[j]] <- net$in.degree[[j]] - 1
  n <- net$outgoing.edges[[i]]
  n <- n[n != j]
  net$outgoing.edges[[i]] <- n
  m <- net$incoming.edges[[j]]
  m <- m[m != i]
  net$incoming.edges[[j]] <- m
  net
}

.ba.remove.node <- function(net, i) {
  n <- net$neighbours[[i]]
  for (j in n) {
    net <- .ba.remove.edge(net, i, j)
  }
  net
}

.ba.add.node.with.edges <- function(net, i, num.edges, offset.exponent, reverse.edges) {
  if (num.edges >= i) stop(sQuote("num.edges"), " must be smaller than ", sQuote("i"))
  
  poss.neighs <- seq_len(i-1) # only smaller than i
  w <- net$degree[poss.neighs]
  w[w == 0] <- 1
  w <- (w / sum(w)) 
  w <- w ^ offset.exponent
  neighs <- sample(poss.neighs, size = num.edges, replace=F, prob=w)
  for (j in neighs) {
    if (reverse.edges) {
      net <- .ba.add.edge(net, j, i)
    } else {
      net <- .ba.add.edge(net, i, j)
    }
  }
  net
}

# adapted from http://stackoverflow.com/questions/24845909/generate-n-random-integers-that-sum-to-m-in-r
#' @importFrom stats rnorm
.ba.num.edges <- function(amnt.nodes, amnt.edges, sd = 1, pos.only = TRUE, redistribute.start = T) {
  vec <- stats::rnorm(amnt.nodes, amnt.edges/amnt.nodes, sd)
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
  if (redistribute.start) {
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
  }
  vec
}

.network.to.df <- function(net) {
  bind_rows(lapply(seq_len(net$amnt.nodes), function(i) {
    bns <- net$outgoing.edges[[i]]
    data.frame(from = rep(i, length(bns)), to = bns)
  }))
}

#' Generate a Barabási–Albert graph
#'
#' @param amnt.nodes number of nodes in the generated graph
#' @param amnt.edges number of edges in the generated graph
#' @param reverse.edges if \code{TRUE} some nodes will have many outgoing links whereas others will not, 
#'   else some nodes will have many incoming links whereas others will not
#' @param offset.exponent the exponent for the weighted sampling. A higher value will result in a 
#'   smaller number of hub nodes having increasingly more outgoing edges.
#' @param trace will generate output in the console if \code{TRUE}
#'
#' @return a list containing several elements:
#' \itemize{
#'   \item amnt.nodes: the number of nodes in the graph
#'   \item degree: the degrees of each node
#'   \item in.degree: the degree of incoming edges of each node
#'   \item out.degree: the degree of outgoing edges of each node
#'   \item incoming.edges: a list of incoming edges for each of the nodes
#'   \item outgoing edges: a list of outgoing edges for each of the nodes
#'   \item data.frame: the edges in the form of a data frame
#' }
#' @export
#' 
#' @importFrom dplyr bind_rows
#'
#' @examples
#' generate.ba(
#'   amnt.nodes = 100,
#'   amnt.edges = 1000, 
#'   reverse.edges = TRUE, 
#'   offset.exponent = 1.5, 
#'   trace = TRUE
#' )
generate.ba <- function(
  amnt.nodes, 
  amnt.edges, 
  reverse.edges = TRUE,
  offset.exponent = 1, 
  trace = FALSE
) {
  if (amnt.edges > amnt.nodes * (amnt.nodes - 1) / 2) 
    stop(sQuote("amnt.edges"), " is too large, as a graph with N nodes can only contain N*(N-1)/2 edges")
  
  net <- .ba.initialise.network(amnt.nodes)
  
  # generate number of edges per node
  nums <- .ba.num.edges(amnt.nodes, amnt.edges)
  
  # add nodes
  for (i in seq_len(amnt.nodes)) {
    if (trace) cat("Generating ", nums[[i]], " ingoing edges for node ", i, "\n", sep="")
    if (nums[[i]] > 0) {
      net <- .ba.add.node.with.edges(net, i, nums[[i]], offset.exponent, reverse.edges = reverse.edges)
    }
  }
  
  # add network as data frame
  df <- bind_rows(lapply(seq_len(net$amnt.nodes), function(i) {
    bns <- net$outgoing.edges[[i]]
    data.frame(i = rep(i, length(bns)), j = bns)
  }))
  net$data.frame <- df
  
  net
}
