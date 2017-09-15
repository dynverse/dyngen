.jacc <- function(seta, setb) {
  aunib <- union(seta, setb)
  if (length(aunib) == 0) {
    0
  } else {
    ain <- aunib %in% seta
    bin <- aunib %in% setb
    sum(ain & bin) / sum(ain | bin)
  }
}
in.common <- function(seta, setb) sum(seta %in% setb)

.ba.mod.add.node.with.edges <- function(net, i, num.edges, exp1, exp2, reverse.edges) {
  #if (num.edges >= i) stop(sQuote("num.edges"), " must be smaller than ", sQuote("i"))
  
  if (num.edges >= 1) {
    poss.neighs <- seq_len(i-1) # only smaller than i
    w <- net$degree[poss.neighs]
    w[w == 0] <- 1
    
    outs <- net$outgoing.edges[poss.neighs]
    
    new.w <- w
    nix <- c()
    pix <- seq_along(poss.neighs)
    while (length(nix) < num.edges) {
      i.sel <- sample(seq_along(pix), size = 1, prob = new.w ^ exp1)
      new.ni <- pix[[i.sel]]
      nix <- c(nix, new.ni)
      if (length(nix) != num.edges) {
        pix <- pix[-i.sel]
        new.w <- new.w[-i.sel]
        zz <- sapply(pix, function(pi) in.common(outs[[pi]], outs[[new.ni]]))
        new.w <- new.w * (log2(zz+1.1) ^ exp2)
      }
    }
    
    neighs <- poss.neighs[nix]
    
    #     i.sel <- sample(seq_along(poss.neighs), size = 1, prob = w ^ exp1)
    #     neighs <- poss.neighs[i.sel]
    #     not.sel.neighs <- poss.neighs[-i.sel]
    #     while (length(neighs) < num.edges) {
    #       if (reverse.edges) {
    #         new.w <- sapply(not.sel.neighs, function(n) sum(sapply(neighs, function(nn) .jacc(net$outgoing.edges[[n]], net$outgoing.edges[[nn]]))))
    #       } else {
    #         new.w <- sapply(not.sel.neighs, function(n) sum(sapply(neighs, function(nn) .jacc(net$outgoing.edges[[n]], net$outgoing.edges[[nn]]))))
    #       }
    #       #new.w[new.w < .01] <- .01
    #       new.w <- new.w + .1
    #       i.sel <- sample(seq_along(not.sel.neighs), size = 1, prob = new.w ^ exp2)
    #       neighs <- c(neighs, not.sel.neighs[[i.sel]])
    #       not.sel.neighs <- not.sel.neighs[-i.sel]
    #     }
    
    for (j in neighs) {
      if (reverse.edges) {
        net <- .ba.add.edge(net, j, i)
      } else {
        net <- .ba.add.edge(net, i, j)
      }
    }
  }
  net
}

#' Generate a Barabási–Albert graph
#'
#' @param amnt.nodes number of nodes in the generated graph
#' @param amnt.edges number of edges in the generated graph
#' @param reverse.edges if \code{TRUE} some nodes will have many outgoing links whereas others will not, 
#'   else some nodes will have many incoming links whereas others will not
#' @param exp1 weight 1
#' @param exp2 weight 2
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
#' generate.ba.with.modules(amnt.nodes = 100, amnt.edges = 1000, reverse.edges = T, exp1 = 1.5, exp2 = 1.5, trace = T)
generate.ba.with.modules <- function(
  amnt.nodes, 
  amnt.edges,
  exp1,
  exp2, 
  reverse.edges = T, 
  trace = F
) {
  if (amnt.edges > amnt.nodes * (amnt.nodes - 1) / 2) 
    stop(sQuote("amnt.edges"), " is too large, as a graph with N nodes can only contain N*(N-1)/2 edges")
  
  net <- .ba.initialise.network(amnt.nodes)
  
  # generate number of edges per node
  nums <- .ba.num.edges(amnt.nodes, amnt.edges, redistribute.start = T)
  
  # add nodes
  for (i in seq_len(amnt.nodes)) {
    if (trace) cat("Generating ", nums[[i]], " ingoing edges for node ", i, "\n", sep="")
    if (nums[[i]] > 0) {
      net <- .ba.mod.add.node.with.edges(net, i, nums[[i]], exp1, exp2, reverse.edges = reverse.edges)
    }
  }
  
  # add network as data frame
  df <- dplyr::bind_rows(lapply(seq_len(net$amnt.nodes), function(i) {
    bns <- net$outgoing.edges[[i]]
    data.frame(i = rep(i, length(bns)), j = bns)
  }))
  net$data.frame <- df
  
  net
}
