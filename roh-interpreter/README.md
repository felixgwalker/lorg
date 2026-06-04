# ROH Interpreter

Interprets runs of homozygosity (ROH) detected in genome-wide SNP or whole-genome sequencing data to infer an individual's or population's inbreeding history, distinguish recent inbreeding from ancient bottleneck effects, and estimate current effective population size. ROH analysis is a cornerstone metric in conservation genomics and de-extinction feasibility assessment, providing a direct measure of genomic inbreeding that integrates both recent and historical demographic signals.

## Approach — Baseline & Novel Layer

**Baseline:** PLINK `--homozyg` and `bcftools roh` detect ROH segments and compute
FROH. This tool does not re-implement those detectors — it calls PLINK/bcftools and
parses their output.

**Novel layer:** PLINK and bcftools assume diploid genotype calls at standard depth.
Ancient DNA genotypes are called pseudo-haploidly (one allele randomly sampled per
site) at coverage < 2×, which systematically inflates false-positive ROH. The novel
contribution is the *aDNA-aware interpretation layer*: a correction model for
pseudo-haploid ROH inflation, a low-coverage confidence filter per ROH segment, and
a de-extinction-specific report framing FROH as a project viability metric (alongside
inbreeding-risk-forecaster output) rather than a standalone statistic.

## Inputs

- Genotype data in VCF format or PLINK binary format (.bed/.bim/.fam)
- Optional: reference population panel in VCF format for comparative context
- Parameters: minimum ROH length thresholds for short (ancient), medium, and long (recent) ROH classes (default: 100 kb, 1 Mb, 10 Mb), minimum SNP density within ROH, sliding window size and step

## Outputs

- A per-individual ROH catalogue in BED format listing all detected ROH intervals with length and SNP count
- An FROH table in CSV format giving genome-wide inbreeding coefficients per individual, stratified by ROH length class
- A ROH length distribution histogram and cumulative FROH barplot (PNG/SVG)
- A demographic inference summary estimating Ne trajectory across historical time periods based on ROH length class proportions
- A population-level summary comparing FROH values and ROH burden across all individuals (CSV and plot)

## Method

Applies a sliding-window homozygosity detection algorithm across phased or unphased genotype data, identifying contiguous runs of homozygous genotype calls above a minimum length threshold and SNP density. Classifies ROH into length bins corresponding to historical time periods: short ROH (>100 kb) reflect ancient bottlenecks (hundreds of generations ago), medium ROH reflect moderate inbreeding (tens of generations), and long ROH (>10 Mb) reflect very recent inbreeding within a few generations. Computes FROH as the sum of ROH length divided by the autosomal genome length. Derives a Ne trajectory by mapping ROH length class proportions onto a coalescent time-scale conversion formula. Optionally benchmarks results against the reference population panel.

## Dependencies

- `pandas` — ROH interval table and FROH summary construction
- `numpy` — sliding-window computations and Ne conversion calculations
- `scipy` — statistical summaries and demographic inference
- `scikit-allel` — VCF and PLINK genotype parsing, ROH detection utilities
- `matplotlib` — ROH length distribution, FROH barplots, and Ne trajectory visualisation
- `pyvcf` — VCF file I/O and filtering
