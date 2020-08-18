# Find possible parent-offspring pairs based on shared alleles

# For each pair of animals, find the loci where they are observed to share at
# least one allele.

# Find all pairs of animal indices
ani_pairs_inds <- combn(n_animals, 2)

# Find indices for first and second animals in each pair
ani_inds_1 <- ani_pairs_inds[1, ]
ani_inds_2 <- ani_pairs_inds[2, ]

# Find pairs with each possible pair of alleles shared at each locus
a1a1_eq_a2a1 <- ales_mat_1[ani_inds_1, ] == ales_mat_1[ani_inds_2, ]
a1a1_eq_a2a2 <- ales_mat_1[ani_inds_1, ] == ales_mat_2[ani_inds_2, ]
a1a2_eq_a2a1 <- ales_mat_2[ani_inds_1, ] == ales_mat_1[ani_inds_2, ]
a1a2_eq_a2a2 <- ales_mat_2[ani_inds_1, ] == ales_mat_2[ani_inds_2, ]

# Find loci where pairs share at least one allele
shared_ales <- a1a1_eq_a2a1 | a1a1_eq_a2a2 | a1a2_eq_a2a1 | a1a2_eq_a2a2



# Find the pairs of animals that share at least one allele at each locus where
# both genotypes are observed.

# This requirement ensures that these pairs are able to be parent and offspring,
# but whether or not that is likely depends on the likely proportion of such
# pairs in the sample, and the information-content of the observed genopairs.
# Other pairs could also be parent and offspring, if genotyping errors or
# mutation have obscured shared alleles.

# Find the missing genopairs at each locus
missing_genopairs <- is.na(shared_ales)

# Find the numbers of loci where both genotypes are observed for each pair
ns_loci_both_obs <- rowSums(!missing_genopairs)

# Find the number of loci where each pair is observed to share at least one allele
ns_loci_ales_shared <- rowSums(shared_ales, na.rm = T)

# Find the indices of possible parent-offspring pairs
possible_pop_inds <- which(ns_loci_ales_shared == ns_loci_both_obs)

# Find the indices of the animals in each parent-offspring pair
possible_pops <- ani_pairs_inds[, possible_pop_inds]



# Save data
save(list = objects(), file = "rdata/whale_analysis_3.Rdata")