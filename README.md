# lorg

A collection of computational biology tools for genomic analysis, de-extinction research, and CRISPR engineering.

*lorg* is Scottish Gaelic for "trace" or "track" — as in following the genomic trail.

## Tools

| Tool | Description |
|------|-------------|
| `adaptive-diversity-scorer` | Scores adaptive genetic diversity in conservation populations by identifying Fst-outlier loci and environment–genotype associations |
| `admixture-signal-scanner` | Scans for admixture signals and estimates per-sample ancestry proportions across a range of K values |
| `allele-frequency-comparator` | Compares allele frequencies across gnomAD populations to identify population-specific variants and quantify differentiation |
| `alternative-splicing-detector` | Detects differential alternative splicing events between two RNA-seq conditions using PSI quantification |
| `ancestral-gene-content-reconstructor` | Reconstructs gene content of ancestral genomes at internal phylogenetic nodes using Dollo parsimony or a Bayesian model |
| `ancestral-state-reconstructor` | Reconstructs ancestral nucleotide or amino acid sequences at internal phylogenetic nodes using marginal ML, parsimony, or Bayesian inference |
| `ancient-dna-damage-classifier` | Classifies DNA damage patterns in ancient/aDNA samples using a Bayesian decay model |
| `ancient-sample-authenticator` | Authenticates ancient DNA samples by evaluating fragment length, deamination, contamination, endogenous fraction, and coverage |
| `annotation-consistency-checker` | Validates a GTF gene annotation for internal consistency: coordinate hierarchy, duplicate IDs, overlapping features, and strand mismatches |
| `assembly-gap-analyser` | Analyses N-run gaps in a genome assembly, classifying their genomic context and potential functional impact |
| `assembly-quality-assessor` | Computes N50, N90, L50, L90, GC content, gap statistics, and ambiguous base count from a genome assembly FASTA |
| `base-edit-outcome-predictor` | Predicts per-base editing probabilities for CBE and ABE base editors using position- and context-dependent efficiency models |
| `biosafety-risk-assessor` | Screens synthetic biology sequences against select agent, virulence factor, and toxin databases to recommend biosafety containment levels |
| `bottleneck-detector` | Detects historical population bottleneck signatures using multiple complementary tests calibrated against coalescent simulations |
| `cas-variant-selector` | Ranks Cas nucleases and editors for a target locus and desired editing goal |
| `chromatin-accessibility-scorer` | Scores chromatin accessibility at genomic peak regions from ATAC-seq BAM data |
| `cnv-significance-assessor` | Scores copy number variants by size, gene content, dosage sensitivity, and population frequency |
| `codon-optimisation-engine` | Optimises protein-coding DNA sequences for expression in a target host using the host codon usage table |
| `coexpression-module-finder` | Identifies coexpression modules from RNA-seq data using WGCNA, clique-based clustering, or k-means |
| `compound-heterozygosity-detector` | Detects compound heterozygous variant pairs in the same gene from phased genotypes or trio data |
| `conservation-priority-ranker` | Ranks populations by conservation urgency using a composite genomic priority score |
| `conserved-synteny-detector` | Detects conserved synteny blocks between two genomes by chaining ortholog anchor pairs |
| `constraint-region-detector` | Detects whether variants fall in genomically constrained regions using gnomAD constraint metrics |
| `contamination-estimator` | Estimates modern human contamination in ancient DNA BAM files using mitochondrial consensus deviation and ANGSD-based approaches |
| `contig-scaffolding-helper` | Orders and orients assembly contigs into scaffolds using paired-end links, Hi-C contacts, or reference-guided scaffolding |
| `crispr-array-designer` | Assembles multi-spacer CRISPR arrays for Cas12a, Cas9, or Cas12b from a list of target sequences |
| `crispr-base-editor-window-visualiser` | Visualises editing windows for base editors across target sequences |
| `crispr-delivery-strategy-selector` | Ranks CRISPR delivery modalities for a given cell type, payload, and experimental context |
| `crispr-knockin-designer` | Designs sgRNA and HDR donor template for CRISPR-mediated precise knock-in |
| `crispr-knockout-designer` | Designs sgRNAs for gene knockout, targeting early coding exons with Doench 2016-style on-target scoring |
| `cross-species-liftover-assistant` | Maps genomic BED intervals from a source species to a target species using a UCSC chain file |
| `demographic-history-inferencer` | Infers demographic history by fitting parametric models to the site frequency spectrum |
| `dna-fragmentation-profiler` | Profiles DNA fragmentation length distribution and post-mortem deamination damage patterns in ancient DNA BAM files |
| `edit-outcome-simulator` | Simulates CRISPR indel outcome distributions at a target site using an inDelphi-style approach |
| `effective-population-size-estimator` | Estimates effective population size from SNP data using Watterson's estimator, LD-based Ne, or a PSMC-style temporal profile |
| `enhancer-conservation-analyser` | Analyses evolutionary conservation of enhancer elements using phastCons/PhyloP scores and multiple alignment coverage |
| `enhancer-target-linker` | Links enhancer elements to candidate target genes using the ABC model, expression correlation, or distance scoring |
| `erv-risk-mapper` | Maps endogenous retroviral element risk in genomic contexts |
| `expression-divergence-scorer` | Scores expression divergence between orthologous gene pairs across two species |
| `founder-effect-estimator` | Estimates the strength and timing of a founder effect by comparing nucleotide diversity and haplotype block length between populations |
| `gene-circuit-stability-estimator` | Estimates dynamic stability and robustness of a synthetic gene circuit by simulating ODE dynamics |
| `gene-essentiality-predictor` | Predicts gene essentiality by aggregating DepMap CRISPR fitness scores, RNAi dependency data, and gnomAD constraint metrics |
| `gene-family-expansion-detector` | Detects gene families that have expanded significantly on specific lineages using a birth-death model |
| `gene-loss-detector` | Detects gene losses on specific lineages by identifying genes absent or pseudogenised relative to outgroups |
| `gene-model-validator` | Validates gene models in a GTF annotation against a reference genome FASTA, checking start/stop codons and splice sites |
| `gene-regulatory-network-builder` | Infers gene regulatory network edges between transcription factors and target genes from expression data |
| `genetic-rescue-candidate-selector` | Ranks donor populations for genetic rescue by balancing heterozygosity gain, kinship distance, and outbreeding depression risk |
| `genetic-rescue-viability-estimator` | Estimates viability outcomes for genetic rescue interventions |
| `genome-completeness-estimator` | Estimates genome assembly completeness by searching BUSCO conserved gene benchmarks against a lineage-specific database |
| `genome-edit-feasibility-scorer` | Scores the overall feasibility of a genome editing project from a composite of locus and context factors |
| `genome-rearrangement-mapper` | Maps chromosomal rearrangements (inversions, translocations, fusions, fissions) between two genomes from synteny block coordinates |
| `genomic-diversity-index` | Computes θW, θπ, Tajima's D, Ho, He, and Fis in sliding windows across the genome |
| `guide-rna-gc-optimiser` | Scores and ranks guide RNAs by GC content features for efficient Cas9-mediated editing |
| `guide-rna-off-target-scorer` | Scores gRNA off-target risk profiles |
| `guide-rna-secondary-structure-analyser` | Analyses secondary structure of sgRNA spacer+scaffold sequences to flag guides with poor accessibility |
| `guide-rna-specificity-ranker` | Ranks sgRNAs by predicted genome-wide specificity using CFD scoring against enumerated off-target sites |
| `hdr-template-designer` | Designs HDR repair templates for precise genome edits |
| `hgt-risk-assessor` | Assesses horizontal gene transfer risk from genomic signals |
| `inbreeding-risk-forecaster` | Forecasts inbreeding risk by estimating current inbreeding from ROH, inferring Ne, and projecting future accumulation |
| `introgression-detector` | Detects gene flow between populations using Patterson's D-statistic and f4-ratio in a four-population ABBA-BABA framework |
| `kinship-coefficient-calculator` | Calculates pairwise kinship coefficients for all samples in a multi-sample VCF and classifies relationships by degree |
| `large-insert-tolerance-predictor` | Predicts tolerance for large genomic insertions |
| `lineage-divergence-dater` | Dates lineage divergence events using a molecular clock, Bayesian dated-tips approach, or pairwise genetic distance |
| `lineage-specific-gene-finder` | Identifies orphan genes in a query species with no detectable homolog in outgroup proteomes |
| `metabolic-pathway-balancer` | Balances a heterologous metabolic pathway for optimal flux distribution using FBA |
| `microhomology-repair-predictor` | Predicts MMEJ deletion products from microhomology sequences flanking a double-strand break |
| `missense-impact-scorer` | Scores the functional impact of missense variants using conservation, substitution cost, and physicochemical property changes |
| `molecular-clock-estimator` | Estimates molecular clock rate and compares strict vs. relaxed clock models using Bayesian MCMC with calibration points |
| `multiplex-edit-planner` | Plans and validates a multiplex CRISPR editing strategy across multiple simultaneous targets |
| `off-target-cluster-detector` | Detects genomic hotspots of CRISPR off-target activity using sliding-window density and clustering |
| `ortholog-mapper` | Maps orthologous genes between a query and target species using reciprocal best BLAST hits or OMA |
| `palaeogenomic-coverage-assessor` | Assesses genome-wide coverage statistics for ancient DNA BAMs: mapping rate, depth, breadth, and duplication rate |
| `pam-flexibility-predictor` | Scores PAM site availability for a panel of Cas variants at a target locus |
| `paralog-cluster-builder` | Builds clusters of paralogous genes from a single-species proteome using all-vs-self BLAST and MCL |
| `peg-rna-optimiser` | Grid-searches PBS and RT template length combinations to find Pareto-optimal pegRNA designs |
| `phylogenetic-distance-estimator` | Estimates phylogenetic distance between species or sequences |
| `population-differentiation-scorer` | Scores population differentiation using Weir-Cockerham Fst, Gst, and Jost's D across all population pairs |
| `population-viability-genomics-estimator` | Estimates population viability by integrating inbreeding, Ne, genetic load, and adaptive diversity to project extinction probability |
| `positive-selection-signal-detector` | Detects positive selection signals in genomic data |
| `post-mortem-damage-simulator` | Simulates post-mortem damage in ancient DNA reads using the Briggs model for deamination and fragment length distributions |
| `prime-edit-design-assistant` | Designs pegRNAs for prime editing from a target locus sequence and desired edit specification |
| `prime-edit-efficiency-predictor` | Predicts the editing efficiency of pegRNA designs using DeepPrime-style sequence features |
| `promoter-strength-estimator` | Estimates promoter strength from sequence features and optional H3K4me3 ChIP-seq signal |
| `promoter-variant-scorer` | Scores variants in promoter regions for transcription factor binding site disruption using position weight matrices |
| `proxy-species-edit-burden-calculator` | Calculates edit burden when using proxy species for de-extinction |
| `rare-variant-prioritiser` | Prioritises rare variants by combining allele frequency, CADD score, gene constraint, and phenotype matching |
| `regulatory-element-conservation-scorer` | Scores conservation of regulatory elements across species |
| `regulatory-rewiring-analyser` | Analyses regulatory rewiring between two species by classifying enhancers and promoters as gained, lost, conserved, or relocated |
| `repair-pathway-bias-estimator` | Estimates NHEJ, MMEJ, and HDR repair pathway probabilities at a CRISPR cut site |
| `repeat-element-classifier` | Classifies repeat elements in a genome assembly by parsing RepeatMasker output |
| `roh-interpreter` | Interprets runs of homozygosity for inbreeding and population analysis |
| `safe-harbour-integration-finder` | Identifies genomic safe harbour sites for stable transgene integration |
| `selection-sweep-detector` | Detects positive selection sweeps using iHS, XP-EHH, CLR, and Tajima's D in sliding windows |
| `splice-impact-predictor` | Predicts the impact of variants on splicing by scoring donor and acceptor sites against position weight matrices |
| `structural-variant-prioritiser` | Prioritises structural variants by breakpoint gene impact, ClinGen dosage sensitivity, and overlap with known pathogenic SVs |
| `synonymous-variant-scorer` | Scores synonymous variants for functional impact |
| `synthetic-promoter-designer` | Designs synthetic promoter sequences by assembling core promoter elements and TFBS arrays |
| `synteny-block-visualiser` | Visualises synteny blocks between genomes |
| `transcription-factor-site-scanner` | Scans FASTA sequences against a JASPAR/MEME PWM database for transcription factor binding sites |
| `utr-variant-analyser` | Analyses variants in 5' and 3' UTR regions for uORF creation/disruption and polyadenylation signal disruptions |
| `variant-pathogenicity-aggregator` | Aggregates evidence from ClinVar and in silico predictors into ACMG/AMP criteria and a five-tier pathogenicity classification |

## Requirements

Each tool has its own `requirements.txt`. Python 3.9+ recommended.

```bash
pip install -r <tool-name>/requirements.txt
```
