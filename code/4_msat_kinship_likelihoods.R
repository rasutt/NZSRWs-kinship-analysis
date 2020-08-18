# Find kinship likelihoods based on microsatellite genotypes

# Find genotype probabilities, and conditional genopair probabilities given that
# the animals are parent and offspring, at each locus

# Create matrix of NAs for genotype probabilities
gt_probs <- matrix(nrow = n_animals, ncol = n_loci)

# Find pairs of animals where second genotype not homozygous
a2a1_not_eq_a2a2 <- ales_mat_1[ani_inds_2, ] != ales_mat_2[ani_inds_2, ]

# Find the number of pairs of animals
n_pairs <- choose(n_animals, 2)

# Create matrix for conditional genotype probabilities.  Need zero because
# adding one component at a time
cond_gt_2_probs_pop <- matrix(0, nrow = n_pairs, ncol = n_loci)

# Loop over loci
for (i in locus_inds) {
  # Find first and second alleles for all animals at this locus 
  ales_mat_1_i <- ales_mat_1[, i]
  ales_mat_2_i <- ales_mat_2[, i]
  
  # Find allele probabilities for first and second alleles for all animals at
  # this locus
  ale_dist_i <- ale_dists[[i]]
  ale_probs_1_i <- ale_dist_i[ales_mat_1_i]
  ale_probs_2_i <- ale_dist_i[ales_mat_2_i]
  
  # Find genotype probabilities for all animals at this locus, the product of
  # the allele probabilities when homozygous, and twice that when heterozygous
  gt_probs[, i] <- 
    ale_probs_1_i * ale_probs_2_i * (1 + (ales_mat_1_i != ales_mat_2_i))
  
  # Find indices of pairs with second genotype first allele shared with first
  # genotype first and second alleles.  which() isn't just faster, it converts
  # NA's to false
  a1a1_eq_a2a1_inds_i <- which(a1a1_eq_a2a1[, i])
  a1a2_eq_a2a1_inds_i <- which(a1a2_eq_a2a1[, i])
  
  # Find indices of pairs with second genotype second allele shared with first
  # genotype first and second alleles and not homozygous
  a1a1_eq_a2a2_a2a1_not_eq_a2a2_inds_i <- which(
    a1a1_eq_a2a2[, i] & a2a1_not_eq_a2a2[, i]
  )
  a1a2_eq_a2a2_a2a1_not_eq_a2a2_inds_i <- which(
    a1a2_eq_a2a2[, i] & a2a1_not_eq_a2a2[, i]
  )
  
  # Find allele probabilities for seconds animals
  ani_2_ale_probs_1_i <- ale_probs_1_i[ani_inds_2]
  ani_2_ale_probs_2_i <- ale_probs_2_i[ani_inds_2]  
  
  # Find conditional probability of second genotype given one allele IBD to one
  # in first genotype
  
  # Each allele in first genotype passed on with probability 0.5.  When allele
  # in second genotype equal add probability of inheriting other by chance
  cond_gt_2_probs_pop[a1a1_eq_a2a1_inds_i, i] <- 
    cond_gt_2_probs_pop[a1a1_eq_a2a1_inds_i, i] +
    0.5 * ani_2_ale_probs_2_i[a1a1_eq_a2a1_inds_i]
  cond_gt_2_probs_pop[a1a2_eq_a2a1_inds_i, i] <- 
    cond_gt_2_probs_pop[a1a2_eq_a2a1_inds_i, i] +
    0.5 * ani_2_ale_probs_2_i[a1a2_eq_a2a1_inds_i]
  
  # Don't add twice when second genotype homozygous
  cond_gt_2_probs_pop[a1a1_eq_a2a2_a2a1_not_eq_a2a2_inds_i, i] <-
    cond_gt_2_probs_pop[a1a1_eq_a2a2_a2a1_not_eq_a2a2_inds_i, i] +
    0.5 * ani_2_ale_probs_1_i[a1a1_eq_a2a2_a2a1_not_eq_a2a2_inds_i]
  cond_gt_2_probs_pop[a1a2_eq_a2a2_a2a1_not_eq_a2a2_inds_i, i] <-
    cond_gt_2_probs_pop[a1a2_eq_a2a2_a2a1_not_eq_a2a2_inds_i, i] +
    0.5 * ani_2_ale_probs_1_i[a1a2_eq_a2a2_a2a1_not_eq_a2a2_inds_i]
}

# Set conditional genopair probabilities for missing data to NA
cond_gt_2_probs_pop[missing_genopairs] <- NA

# Remove pairs with second genotype homozygous, and missing genopairs, as large
# and no longer needed
rm(a2a1_not_eq_a2a2, missing_genopairs)



# Find genopair probabilities given that the animals are urelated, parent and
# offspring, and self, at each locus

# Find genopair probabilities at all loci given an unrelated pair
gp_probs_up <- gt_probs[ani_inds_1, ] * gt_probs[ani_inds_2, ]

# Find genopair probabilities given a parent-offspring pair, the product of the
# probability of the first genotype and the conditional probability of the
# second genotype, given that one allele is IBD
gp_probs_pop <- gt_probs[ani_inds_1, ] * cond_gt_2_probs_pop

# Remove conditional probabilities of second genotype given a parent-offspring
# pair, as large and no longer needed
rm(cond_gt_2_probs_pop)

# Find loci where both alleles shared, and their indices, for all pairs of
# animals
both_ales_shared <- 
  (a1a1_eq_a2a1 & a1a2_eq_a2a2) | (a1a1_eq_a2a2 & a1a2_eq_a2a1) 
both_ales_shared_inds <- which(both_ales_shared)

# Create matrix of NAs for genopair probabilities given a self-pair
gp_probs_sp <- matrix(nrow = n_pairs, ncol = n_loci)

# Find genopair probabilities given a self-pair
gp_probs_sp[both_ales_shared_inds] <- 
  gt_probs[ani_inds_1, ][both_ales_shared_inds]
gp_probs_sp[which(!both_ales_shared)] <- 0

# Remove pairs with different pairs of alleles shared at each locus, and loci
# where both alleles shared, as large and no longer needed
rm(a1a1_eq_a2a1, a1a2_eq_a2a2, a1a1_eq_a2a2, a1a2_eq_a2a1
   # both_ales_shared
   )



# Find log genopair probabilities given simple kinships

# Find log genopair probabilities given an unrelated pair
kin_log_likes_msats <- data.frame(UP = rowSums(log(gp_probs_up), na.rm = T))

# Find log genopair probabilities given a half-sibling pair
kin_log_likes_msats$SOK <- rowSums(
  log(1/2 * gp_probs_up + 1/2 * gp_probs_pop),
  na.rm = T
)

# Find log genopair probabilities given a parent-offspring pair
kin_log_likes_msats$OP <- kin_log_likes_msats$PO <- 
  kin_log_likes_msats$PO_OP <- rowSums(log(gp_probs_pop), na.rm = T)

# Find log genopair probabilities given a full-sibling pair
kin_log_likes_msats$SOKx2 <- rowSums(
  log(1/4 * gp_probs_up + 1/2 * gp_probs_pop + 1/4 * gp_probs_sp),
  na.rm = T
)

# Find log genopair probabilities given a self-pair
kin_log_likes_msats$SP <- rowSums(log(gp_probs_sp), na.rm = T)

# Remove genopair probabilities given kinships, as large and no longer needed
# for a while
rm(gp_probs_up, gp_probs_pop, gp_probs_sp)



# Save data as R objects
save(list = objects(), file = "rdata/whale_analysis_4.Rdata")