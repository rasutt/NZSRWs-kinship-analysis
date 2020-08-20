# Function to find probabilities for all possible genotypes at one locus

# The currect version of the code finds each genotype probability for each
# animal individually, rather than finding all possible, as below, and then
# indexing for each animal.  I was thinking it was slow and complicated because
# of the need to flatten it to index it in a vectorised way, but now I know I
# can index by coordinates, so it could end up being nice.
find_genotype_probs_one_locus <- function(ale_dists_one_locus) {
  # Find products of allele probabilities for all possible genotypes 
  outer_mat <- outer(ale_dists_one_locus, ale_dists_one_locus)
  
  # Double them for heterozygous genotypes and index just one set?
  full_lower_tri_inds <- lower.tri(outer_mat, diag = T)
  genotype_probs_one_locus <- 
    ((full_lower_tri_inds + lower.tri(outer_mat)) *
       outer_mat)[full_lower_tri_inds]
}
