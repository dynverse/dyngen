# dyngen: Synthetic single cell data

Simulate different types of experiments:

* Cellular heterogeneity
* Normal trajectories
* Branched trajectories
* Cycles
* Multiple parallel dynamic processes and their interactions
* Comparing trajectories between
  * Patients
  * Genotypes
  * Similar cell types in different environments
* Integrating different data types prior or after trajectory inference
* Synchronizing trajectories between interacting cells (physically and/or chemically)
* Dynamics of the regulatory network

Get data both for individual cells as well as at the bulk level

Get different data types:

* Gene expression
* Protein expression

Stress test your methods against:

* Noise
* Contamination

## Notes

* Do we need to handle transcriptional burts? (eg. see last paragraph of online methods in http://www.nature.com/nmeth/journal/v11/n6/full/nmeth.2930.html). This would dramatically increase the stochasticity of gene expression.
* Do we need to model individual binding events? I think we should, but this would greatly increase the number of reactions... (and not sure whether it will really have a big effect)

## Interesting studies

* http://genome.cshlp.org/content/24/3/496.full