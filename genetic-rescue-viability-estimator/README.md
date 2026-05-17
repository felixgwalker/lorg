# Genetic Rescue Viability Estimator

Estimates the likely success and population-genetic outcome of a genetic rescue intervention in an endangered or functionally extinct population, modelling how introduction of new genetic material will reduce inbreeding depression, restore heterozygosity, and improve long-term viability. This tool helps conservation managers and de-extinction teams evaluate whether genetic rescue is a viable near-term intervention and model expected outcomes before committing resources to translocation or assisted reproduction programmes.

## Inputs

- Population genotype data in VCF format (all sampled individuals)
- Optional pedigree file or kinship matrix in CSV format
- Demographic parameters: current population size (N), target population size, number of rescue individuals to introduce, number of generations to simulate
- Optional: donor population VCF if distinct from the recipient population
- Parameters: selection coefficient against inbred individuals, number of simulation replicates

## Outputs

- Inbreeding coefficient (F) estimates per individual in CSV format
- Projected heterozygosity and inbreeding trajectories over the specified generation range (CSV and plot)
- A genetic rescue viability score summarising the expected magnitude of improvement
- A visualisation of pre- and post-rescue heterozygosity distributions and inbreeding decline (PNG/SVG)
- A simulation summary report with confidence intervals across replicates

## Method

Computes individual inbreeding coefficients using both ROH-based (FROH) and genomic relatedness matrix approaches. Estimates current heterozygosity and effective population size (Ne). Simulates forward-in-time population genetics under a Wright-Fisher model incorporating the rescue individuals' genotypes. Models inbreeding depression as a fitness function of F. Tracks heterozygosity, F, and population persistence probability across generations. Aggregates replicates to produce mean trajectories and confidence intervals for the viability score.

## Dependencies

- `pandas` — genotype table handling and output formatting
- `numpy` — Wright-Fisher simulation and matrix operations
- `scipy` — statistical summaries and confidence interval computation
- `scikit-allel` — VCF parsing, ROH detection, and population genetic statistics
- `matplotlib` — trajectory and distribution visualisations
- `networkx` — pedigree graph construction and relatedness computation
