# Run analysis from any step.  Results are saved at the end of each script

source("code/1_import_and_process_msats.R")
# load(file = "rdata/whale_analysis_1.Rdata")
source("code/2_add_covariates_and_haplotypes.R")
# load(file = "rdata/whale_analysis_2.Rdata")
source("code/3_possible_pops_from_shared_ales.R")
# load(file = "rdata/whale_analysis_3.Rdata")
source("code/4_msat_kinship_likelihoods.R")
# load(file = "rdata/whale_analysis_4.Rdata")
source("code/5_msat_haplotype_sex_kinship_likelihoods.R")
# load(file = "rdata/whale_analysis_5.Rdata")
source("code/6_kinship_analyses.R")
# load(file = "rdata/whale_analysis_6.Rdata")
