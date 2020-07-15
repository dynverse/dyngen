library(tidyverse)
library(dyngen)

set.seed(1)

### Step 1: Define backbone and other parameters
# A dyngen simulation can be started by providing a backbone to the `initialise_model()` function.
# The backbone of a `dyngen` model is what determines the overall dynamic process 
# that a cell will undergo during a simulation. It consists of a set of gene modules, which regulate
# eachother in such a way that expression of certain genes change over time in a specific manner.
model <- 
  initialise_model(
    num_tfs = 12,
    num_targets = 30,
    num_hks = 15,
    backbone = backbone_bifurcating(),
    verbose = TRUE,
    simulation_params = simulation_default(
      census_interval = 10,
      ssa_algorithm = ssa_etl(tau = 300 / 3600)
    ),
    download_cache_dir = "~/.cache/dyngen",
    # be sure to change this to the number of cores on your computer
    # when in doubt, set to a value of 4
    num_cores = 8 
  )

plot_backbone_statenet(model)
plot_backbone_modulenet(model)

# For backbones with all different sorts of topologies, check `list_backbones()`:
names(list_backbones())

### Step 2: Generate transcription factors (TFs)
# Each gene module consists of a set of transcription factors.
# These can be generated and visualised as follows.
model <- generate_tf_network(model)
plot_feature_network(model, show_targets = FALSE)

### Step 3: Sample target genes and housekeeping genes (HKs) 
# Next, target genes and housekeeping genes are added to the network by
# sampling a gold standard gene regulatory network using the Page Rank algorithm.
# Target genes are regulated by TFs or other target genes, while HKs are only regulated
# by themselves.
model <- generate_feature_network(model)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)

### Step 4: Generate kinetics
# Note that the target network does not show the effect of some interactions, 
# because these are generated along with other kinetics parameters of the 
# SSA simulation.
model <- generate_kinetics(model)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)

### Step 5: Simulate gold standard
# The gold standard is simulated by enabling certain parts of 
# the module network and performing ODE simulations. The gold standard
# are visualised by performing a dimensionality reduction on the 
# mRNA expression values.
model <- generate_gold_standard(model)
plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2")

# The expression of the modules (average of TFs) can be visualised as follows.
plot_gold_expression(model, what = "mol_mrna") # mrna
plot_gold_expression(model, label_changing = FALSE) # premrna, mrna, and protein

### Step 6: Simulate cells.
# Cells are simulated by running SSA simulations. The simulations are again
# using dimensionality reduction.
model <- generate_cells(model)
plot_simulations(model)

# The gold standard can be overlayed on top of the simulations.
plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2")

# We can check how each segment of a simulation is mapped to the gold standard.
plot_gold_mappings(model, do_facet = FALSE) + scale_colour_brewer(palette = "Dark2")

# The expression of the modules (average of TFs) of a single simulation can be visualised as follows.
plot_simulation_expression(model, 1:4, what = "mol_mrna")

### Step 7: Experiment emulation
# Effects from performing a single-cell RNA-seq experiment can be emulated as follows.
model <- generate_experiment(model)

### Step 8: Convert to a dynwrap object (optional)
dataset <- wrap_dataset(model)

### Visualise with `dyno`
library(dyno)
plot_dimred(dataset)
plot_graph(dataset)

### Infer trajectory on expression data
library(dyno)
pred <- infer_trajectory(dataset, ti_slingshot())
plot_dimred(pred)




## One-shot function
# `dyngen` also provides a one-shot function for running 
# all of the steps all at once and producing plots.
set.seed(1)

model <- 
  initialise_model(
    num_tfs = 12,
    num_targets = 30,
    num_hks = 15,
    backbone = backbone_bifurcating_converging(),
    verbose = TRUE,
    simulation_params = simulation_default(
      census_interval = 10,
      ssa_algorithm = ssa_etl(tau = 300 / 3600)
    ),
    download_cache_dir = "~/.cache/dyngen",
    # be sure to change this to the number of cores on your computer
    # when in doubt, set to a value of 4
    num_cores = 8 
  )

out <- generate_dataset(
  model,
  make_plots = TRUE
)

dataset <- out$dataset
model <- out$model

# write output to files for later reuse
write_rds(dataset, "dataset.rds", compress = "gz")
write_rds(model, "model.rds", compress = "gz")
ggsave("plot.pdf", out$plot, width = 20, height = 15)

# `dataset` and `model` can be used in much the same way as before.
plot_dimred(dataset)
plot_graph(dataset)
pred <- infer_trajectory(dataset, ti_slingshot(), verbose = FALSE)
plot_dimred(pred)

