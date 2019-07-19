# DNA-based-predictors-of-non-genetic-cancer-phenotypes

Analysis scripts and relative data used in the paper "Genetic determinants of the molecular portraits of epithelial cancers" by Xia et al. 

## Data
contains many .rda files ready to use in the analysis.

## Rscripts
**helper.R**: This script has many R functions used in the analysis.

**signature_score_and_segment_score_calculation.R**: This script is used to calculate signature scores from RNA expression data and segment scores from gene-level copy number data.

**association_test.R**: This script is used to perform genome wide association tests between gene signatures and DNA copy number alterations. 

**Elastic_Net_modeling.R**: This script is used to perfrom Elastic Net modeling analysis. All data used to build Elastic Net models are included in Data folder. 