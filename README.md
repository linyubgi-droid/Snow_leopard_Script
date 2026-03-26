# Snow_leopard_Script

Scripts for genomic insights on conservation priorities for global snow leopards.

## Directory structure and description

### 01_Mapping
Read alignment and preprocessing of raw sequencing data, including quality control, mapping to the reference genome, sorting, and duplicate removal.

### 02_Joint_call
Variant calling pipeline, including individual GVCF generation and joint genotyping to produce a unified variant dataset across all samples.

### 03_PCA
Principal Component Analysis (PCA) to infer population structure based on genome-wide SNP data.

### 04_TREE
Phylogenetic tree construction using genome-wide variants to infer evolutionary relationships among individuals or populations.

### 05_Admixture
Population structure analysis using ADMIXTURE, including estimation of ancestry proportions and cross-validation to determine optimal K.

### 06_HET
Genome-wide heterozygosity estimation for each individual.

### 07_ROH
Detection of runs of homozygosity (ROH) to assess inbreeding levels and demographic history.

### 08_mutation_load
Estimation of mutational load, including annotation and classification of deleterious variants.

### 09_treemix
Inference of population splits and gene flow using TreeMix.

### 10_Dsuite
D-statistics (ABBA-BABA related tests) using Dsuite to detect introgression signals.

### 11_ABBA-BABA
Additional ABBA-BABA analyses to quantify gene flow among populations.

### 12_MSMC2
Demographic history inference using MSMC2, including effective population size changes over time.

### 13_PopSizeABC
Approximate Bayesian Computation (ABC)-based estimation of population size history.

### 14_qpgraph
Admixture graph modeling using qpGraph to test complex evolutionary scenarios.

### 15_fastsimcoal
Demographic modeling using fastsimcoal to infer population parameters under different evolutionary models.

### 16_SLIM
Forward simulations using SLiM to model evolutionary processes and validate demographic scenarios.

### 17_Climate_adaption
Analysis of climate adaptation signals, including environmental association and genomic offset estimation.

---

## Notes

All scripts are designed to reproduce the analyses described in the manuscript. Users should adapt file paths and parameters according to their own datasets.
