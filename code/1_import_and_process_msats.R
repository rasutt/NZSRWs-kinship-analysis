# Code to import and process NZ southern right whale microsatellite genotype
# data.

# Read in data from single-sheet file containing one sample-name and the
# microsatellite genotype for each animal sampled.
whale_data <- read.csv(
  "spreadsheet_csvs/All_NZ_CI_AI_MNZ_incl_calves_cervus_with_2019genotypes - Sheet1.csv"
)

# Read data from 2009 sheet from big multi-sheet file
sheet_2009 <- read.csv(
  "spreadsheet_csvs/SRW_database_V2_July2012_NoCalves_repaired.xlsx - 2009.csv"
)

# Copy missing data from 2009 sheet for animals in table of all unique animals.
# Looks like the microsatellite genotypes got shifted to the left seven places
# and cut off
whale_data[
  whale_data$Sample.Name %in% c("Eau09AI162", "Eau09AI208", "Eau09AI215"),
  2:27
  ] <- sheet_2009[
  sheet_2009$Sample.Name %in% c("Eau09AI162", "Eau09AI208", "Eau09AI215"),
  6:31
  ]



# Find metadata and process data for easy access

# Find the numbers of animals and loci
n_animals <- nrow(whale_data)
n_loci <- (ncol(whale_data) - 1) / 2

# Find locus indices
locus_inds <- 1:n_loci

# Find first and second allele indices
ale_inds_1 <- 2 * locus_inds - 1
ale_inds_2 <- 2 * locus_inds

# Find locus names
ale_names_1 <- names(whale_data)[ale_inds_1 + 1]
locus_names <- substr(ale_names_1, 1, nchar(ale_names_1) - 1)

# Convert data to matrix
ales_mat <- as.matrix(whale_data[, -1])

# Find first and second alleles for each animal at each locus
ales_mat_1 <- ales_mat[, ale_inds_1]
ales_mat_2 <- ales_mat[, ale_inds_2]



# Find allele distributions and convert allele data to factors with matching indices

# Create list for allele distributions
ale_dists <- vector("list", n_loci)

# Loop over the loci
for (i in locus_inds) {
  # Find distribution of observed alleles at this locus
  ale_tab <- table(
    c(ales_mat_1[, i], ales_mat_2[, i]),
    exclude = c(NA, "0", "  0", "OOB")
  )
  ale_dists[[i]] <- ale_tab / sum(ale_tab)
  
  # Convert allele columns to factors with indices matching distribution over
  # both columns.   Missing data is converted to NA
  whale_data[, ale_inds_1[i] + 1] <- factor(
    ales_mat_1[, i],
    levels = names(ale_dists[[i]])
  )
  whale_data[, ale_inds_2[i] + 1] <- factor(
    ales_mat_2[, i],
    levels = names(ale_dists[[i]])
  )
}



# Save data as R objects
save(list = objects(), file = "rdata/whale_analysis_1.Rdata")