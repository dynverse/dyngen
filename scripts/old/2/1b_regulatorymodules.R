modulenet_to_modules = function(modulenet, modulenodes, ngenespermodule=4) {
  ngenespermodule = 4
  modulemembership = lapply(c(modulenet$from, modulenet$to) %>% unique, function(x) ((x-1)*ngenespermodule):((x-1)*ngenespermodule+ngenespermodule-1))
  net = lapply(1:nrow(modulenet), function(i) {
    from_module = modulenet[i,]$from
    to_module = modulenet[i,]$to
    from = modulemembership[[from_module]]
    to = modulemembership[[to_module]]
    nedges = length(to) * 0.5
    
    effect = modulenet[i,]$effect
    strength = modulenet[i, ]$strength
    
    net = expand.grid(from=from, to=to) %>% as_tibble() %>% group_by(to) %>% sample_n(nedges) %>% mutate(effect=effect, strength=strength, cooperativity=2)
    #tibble(from=sample(from, nedges, T), to=sample(to, nedges, T), effect=modulenet[i,]$effect, strength=modulenet[i, ]$strength, cooperativity=2)
  }) %>% bind_rows
  
  tibble::lst(modulemembership, net)
}

## add extra target genes to every ldtf
add_targets_individually = function(net) {
  allgenes = ldtfs = sort(unique(union(net$from, net$to)))
  for(ldtf in ldtfs) {
    nnewtargets = sample(1:2, 1)
    if (nnewtargets > 10) {
      subnet = dyngen::generate.ba.with.modules(nnewtargets, nnewtargets*2, 0.05, 0.05)$data.frame %>% rename(from=i, to=j) %>% mutate(effect=0, strength=1, cooperativity=1) %>% mutate(from=from+max(allgenes), to=to+max(allgenes)) %>% mutate(from=replace(from, from==max(allgenes)+1, ldtf))
    } else if (nnewtargets > 0) {
      subnet = tibble(from=ldtf, to=(max(allgenes)+1):(max(allgenes) + nnewtargets + 1), effect=0, strength=1, cooperativity=1)
    } else {
      subnet= tibble()
    }
    
    net = bind_rows(net, subnet)
    allgenes = sort(unique(union(net$from, net$to)))
  }
  
  tibble::lst(net, allgenes)
}

## add extra target genes to every ldtf in a modular nature, where some target modules are regulated by multiple ldtf modules
add_targets_shared = function(net) {
  allgenes = ldtfs = sort(unique(union(net$from, net$to)))
  for(i in 1:length(ldtfs)) {
    nnewtargets = sample(1:6, 1)
    newtargets = (max(allgenes)+1):(max(allgenes) + nnewtargets + 1)
    
    subnet = expand.grid(from=sample(ldtfs, size = sample(1:3, size=1, prob=c(0.4, 0.4, 0.2))), to=newtargets) %>% as.data.frame()
    
    subnet = subnet %>% mutate(effect=0, strength=1, cooperativity=1)
    
    net = bind_rows(net, subnet)
    allgenes = sort(unique(union(net$from, net$to)))
  }
  tibble::lst(net, allgenes)
}

