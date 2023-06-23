# ENCODE4 distal regulation CRISPR data
Re-analyze and harmonize Gasperini et al., 2019, Schraivogel et al., 2020 and Fulco et al., 2019 CRISPRi enhancer screens for ENCODE distal regulation E-G benchmarking.

Re-analyzes scRNA-seq based screens (Gasperini et al., 2019 and Schraivogel et al., 2020) by performing differential expression tests to identify enhancer-gene pairs using MAST. Also applies a simulation-based power analysis to estimate the statistical power of each tested enhancer-gene pair to detect certain perturbation effect sizes (e.g. 80% power to detect a 15% decrease in expression upon perturbation).

Filters re-analyzed data based on statistical power and performs a liftover from hg19 to GRCh38 for Gasperini et al., 2019. Combines re-analyzed data with downloaded Fulco et al., 2019 data and creates output files containing element-gene pairs used for benchmarking predictive models by the ENCODE4 NAWG distal regulation subgroup.
