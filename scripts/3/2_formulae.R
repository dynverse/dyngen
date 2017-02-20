generate_formulae = function(net) {
  ## generating reaction formulas (the probability that a reaction occurs in [t, t+dt])
  formulae = list()
  nus.changes = list()
  add.formula = function(formula, nu) {formulae <<- c(formulae, formula);nus.changes <<- c(nus.changes, list(nu))}
  variables = list() # store every variable including some extra info for later use
  vargroups = list() # store every variable in groups, for example if we want all k values later on
  add.variable = function(variable, ...) {
    varname = variable@string
    info = list(name=varname, ...)
    variables[[varname]] <<- info
    
    group = str_replace(varname, "^(.*?)[_$].*", "\\1")
    if (group %in% c("b", "u", "x", "y", "k", "r", "d", "p", "q", "a0", "a")) {
      vargroups[[group]] <<- c(vargroups[[group]], varname)
    }
    return(variable)
  }
  get.nu = function(higher=NULL, lower=NULL) {setNames(c(rep(1, length(higher)), rep(-1, length(lower))), lapply(c(higher, lower), function(x) {x@string}))}
  
  kg = add.variable(fvar("kg")) # global tf binding speed (influences burstiness of transcriptional regulation)
  rg = add.variable(fvar("rg")) # global transcription rate
  dg = add.variable(fvar("dg")) # global RNA degradation rate
  pg = add.variable(fvar("pg")) # global translation rate
  qg = add.variable(fvar("qg")) # global protein degradation rate
  for (target in G) {
    x = add.variable(fvar("x", target), gene=target)
    y = add.variable(fvar("y", target), gene=target)
    r = add.variable(fvar("r", target), gene=target)
    d = add.variable(fvar("d", target), gene=target)
    q = add.variable(fvar("q", target), gene=target)
    p = add.variable(fvar("p", target), gene=target)
    a0 = add.variable(fvar("a0", target), gene=target)
    
    subnet = net[net$to == target,]
    regs = subnet$from
    effects = subnet$effect
    
    geneinfo = genes[genes$gene == target,] %>% as.list
    
    effects[effects == 0] = sample(c(1, -1), sum(effects == 0), replace=T)
    net[net$to == target,]$effect = effects
    
    boundvars = list()
    
    if(!is.na(geneinfo$a0)) {
      a0_value = geneinfo$a0
    } else if(length(regs) == 0) {
      a0_value = 1
    } else if(sum(effects == -1) == 0) {
      a0_value = 0.0001
    } else if(sum(effects == 1) == 0) {
      a0_value = 1
    } else {
      a0_value = 0.5
    }
    variables[[a0@string]]$value = a0_value
    
    # regulator binding sites
    if (length(regs) > 0) {
      for (i in seq_len(length(regs))) {
        regulator = regs[[i]]
        effect = effects[[i]]
        strength = subnet[i,]$strength
        cooperativity = subnet[i,]$cooperativity
        
        #u = add.variable(fvar("u", target, regulator), target=target, regulator=regulator)
        #b = add.variable(fvar("b", target, regulator), target=target, regulator=regulator)
        k = add.variable(fvar("k", target, regulator), target=target, regulator=regulator, strength=strength)
        
        #boundvars = c(boundvars, b)
        
        y_regulator = fvar("y", regulator)
        
        # binding
        if (cooperativity > 1) {
          #formula = fcon(paste0(u@string, " * ", y_regulator@string, "^", cooperativity, " * ", kg@string))
          #formula = u * y_regulator * y_regulator * kg
        } else {
          #formula = u * y_regulator * kg
        }
        
        #nu = get.nu(b, list(y_regulator, u))
        #add.formula(formula, nu)
        
        # deassociation
        #formula = b * k * kg
        #nu = get.nu(list(y_regulator, u), b)
        #add.formula(formula, nu)
        
        # add to decision tree
        a = add.variable(fvar("a", target, regulator), target=target, regulator=regulator, effect=effect)
        
        # quick and dirty decision tree, the last interaction will have priority
        if(i == 1) {
          decisiontree = fcon(paste0(a0@string, "+", a@string, " * (", y_regulator@string, "/", k@string, ")^", cooperativity))
          decisiontree_denom = fcon(paste0(1, "+ (", y_regulator@string, "/", k@string, ")^", cooperativity))
          
          #decisiontree = fcon(paste0(a0@string, "+", b@string, "/", u@string, "*", a@string))
          #decisiontree_denom = fcon(paste0(1, "+", b@string, "/", u@string))
          
          #decisiontree = fcon(paste0(a0@string, "+", b@string, "*", a@string))
          #decisiontree_denom = fcon(paste0(1, "+", b@string))
        } else {
          decisiontree = fcon(paste0(decisiontree@string, "+", a@string, " * (", y_regulator@string, "/", k@string, ")^", cooperativity))
          decisiontree_denom = fcon(paste0(decisiontree_denom@string, "+ (", y_regulator@string, "/", k@string, ")^", cooperativity))
          
          #decisiontree = fcon(paste0(decisiontree@string, "+", b@string, "/", u@string, "*", a@string))
          #decisiontree_denom = fcon(paste0(decisiontree_denom@string, "+", b@string, "/", u@string))
          
          #decisiontree = fcon(paste0(decisiontree@string, "+", b@string, "*", a@string))
          #decisiontree_denom = fcon(paste0(decisiontree_denom@string, "+", b@string))
        }
      }
    } else {
      decisiontree = a0
      decisiontree_denom = fcon(1)
    }
    
    decisiontree = fcon(paste0("(", decisiontree@string, ")/(", decisiontree_denom@string, ")"))
    
    # mRNA production
    formula = rg * r * decisiontree
    nu = get.nu(x)
    add.formula(formula, nu)
    
    # mRNA degradation
    formula = dg * d * x
    nu = get.nu(lower=x)
    add.formula(formula, nu)
    
    if(target %in% net$from) {
      # protein production
      formula = pg * p * x
      nu = get.nu(y)
      add.formula(formula, nu)
      
      # protein degradation
      formula = qg * q * y
      nu = get.nu(lower=y)
      add.formula(formula, nu)
    }
  }
  
  
  ###
  formula = fcon(1)
  nu = get.nu(higher=rg)
  #add.formula(formula, nu)
  
  formula = fcon(1)
  nu = get.nu(lower=rg)
  #add.formula(formula, nu)
  
  formula = (x / 10)
  nu = get.nu(lower=rg)
  #add.formula(formula, nu)
  
  #"
  
  formulae.strings <- sapply(formulae, function(fl) fl@string)
  
  named.list(vargroups, variables, formulae.strings, nus.changes)
}