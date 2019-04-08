names <- list.files("inst/ext_data/modulenetworks/")

quote <- function(x) paste0("\"", x, "\"")

strs <- map_chr(names, function(mod) {
  cat(mod, "\n", sep = "")
  edge_operations <- read_tsv(glue::glue("inst/ext_data/modulenetworks/{mod}/edge_operations.tsv"))
  module_nodes <- read_tsv(glue::glue("inst/ext_data/modulenetworks/{mod}/modulenodes.tsv"))
  module_net <- read_tsv(glue::glue("inst/ext_data/modulenetworks/{mod}/modulenet.tsv"))
  
  paste0(
    "#' @export\n",
    "#' @rdname module_network\n",
    "module_network_", mod, " <- function() {\n",
    "  nodes <- tribble(\n",
    "    ~module_id, ~a0, ~burn,\n",
    "    ",
    module_nodes %>%
      mutate_at(vars(module_id), quote) %>% 
      select(module_id, a0, burn) %>% 
      mapdf_chr(paste, collapse = ", ") %>% 
      paste(collapse = ",\n    "),
    "\n",
    "  )\n",
    "\n",
    "  edges <- tribble(\n",
    "    ~from, ~to, ~effect, ~strength, ~cooperativity,\n",
    "    ",
    module_net %>%
      mutate_at(vars(from, to), quote) %>% 
      select(from, to, effect, strength, cooperativity) %>% 
      mapdf_chr(paste, collapse = ", ") %>% 
      paste(collapse = ",\n    "),
    "\n",
    "  )\n",
    "\n",
    "  operations <- tribble(\n",
    "    ~from, ~to, ~module_progression, ~start, ~burn,\n",
    "    ",
    edge_operations %>% 
      mutate_at(vars(from, to, module_progression), quote) %>% 
      select(from, to, module_progression, start, burn) %>% 
      mapdf_chr(paste, collapse = ", ") %>% 
      paste(collapse = ",\n    "),
    "\n",
    "  )\n",
    "  \n",
    "  module_network(nodes, edges, operations)\n",
    "}\n",
    "\n"
  )
})

write_lines(unlist(strs), "R/1a_module_networks_gen.R")
