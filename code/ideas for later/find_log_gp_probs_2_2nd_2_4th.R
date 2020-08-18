# Find genopair probabilities given two second-order, and two fourth-order
# kinships.  Full-sibs, and parents half-sibs, so also twice first-cousins.

# The following is wrong because it doesn't take into account homozygosity
# within each animal, which is a factor when parents are related.  But it might
# be part of the right expression?

# Add the probabilities of the allele you get from your mum being the one they
# get from their dad, because they're the same one from your grandmum, and vice
# versa.  Each has probability 1/8, giving (1/8)^2 for both at once, and 2 * 1/8
# * 7/8 for one at a time.  Then just remember not to double count for the case
# of one allele ibd.
log_gp_probs_2_2nd_2_4th <- rowSums(
  log(
    (1/4 - 15/64 + 2 * 1/4 * 7/64) * gp_probs_up + 
      (1/2 + 14/64 - 2 * 1/4 * 7/64) * gp_probs_pop + 
      (1/4 + 1/64) * gp_probs_sp
  )
)

# Combine kinship likelihoods for likely parent-offspring pairs according to
# shared alleles
kinship_likelihoods <-  
  cbind(
    log_gp_probs_up,
    log_gp_probs_hsp,
    log_gp_probs_pop, 
    log_gp_probs_fsp, 
    log_gp_probs_sp,
    log_gp_probs_2_2nd_2_4th
  )[is.finite(log_gp_probs_pop), ]

# Show first rows
head(kinship_likelihoods)

# Find number of likely parent-offspring pairs according to shared alleles that
# are more likely to be related by two second-order and two fourth-order
# kinships
sum(kinship_likelihoods[3, ] > kinship_likelihoods[1, ])