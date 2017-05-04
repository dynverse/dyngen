generate_random_tree = function() {
  dig_deeper = function(net=tibble(from=c(1), to=c(2)), level=0, parentstate=2, decay=1.4) {
    if(decay <= 1) stop("decay has to be > 1")
    for (i in c(1,2)) {
      maxstate = ifelse(nrow(net), max(net$to), 0)
      if(runif(1) < decay^(-level)) {
        net = net %>% bind_rows(tibble(from=parentstate, to=maxstate+1))
        net = dig_deeper(net, level+1, maxstate+1, decay)
      }
    }
    return(net)
  }
  
  states = dig_deeper()
  states = states %>% group_by(from) %>% summarize(singular=n() == 1) %>% right_join(states, by="from")
  
  states = 
  
  igraph::graph_from_data_frame(states[, c("from", "to")]) %>% plot
}
