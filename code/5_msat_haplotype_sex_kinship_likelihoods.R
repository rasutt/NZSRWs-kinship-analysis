# Find kinship likelihoods based on MtDNA haplotypes and sexes, and combine with
# those based on microsatellite genotypes

# The former are sometimes ignored when the haplotypes are known but one or both
# of the sexes are unknown (eg. when it is less trivial than for unrelated
# pairs) but this almost never happens in this data. It would be easy to add
# later, as they are just the weight averages of the possible cases.

# The sex probabilities are approximated by 0.5, although it is more like 0.55
# for being female.  But these are only relevant to distinguishing self-pairs,
# which is easy with any decent genotype

# Prepare data for finding conditional probabilities for haplotypes given
# kinships and sexes

# Find haplotype distribution
ht_tab <- table(whale_data$DLP)
ht_dist <- ht_tab / sum(ht_tab)

# Find haplotype probabilities for first and second animals
ht_probs <- ht_dist[whale_data$DLP]
ht_probs_1 <- ht_probs[ani_inds_1]
ht_probs_2 <- ht_probs[ani_inds_2]

# Find sex combinations for all pairs.  Only M and F should be allowed
same_sex <- whale_data$SEX[ani_inds_1] == whale_data$SEX[ani_inds_2]
both_fem <- same_sex & (whale_data$SEX[ani_inds_1] == "F")
both_male <- same_sex & !both_fem
male_1_fem_2 <- !same_sex & (whale_data$SEX[ani_inds_1] == "M")
fem_1_male_2 <- !same_sex & !male_1_fem_2

# Find whether haplotypes equal for all pairs
same_hts <- whale_data$DLP[ani_inds_1] == whale_data$DLP[ani_inds_2]

# Find indices of combinations of sexes and equal haplotypes for all pairs where
# both known.  which() is both faster and treats NA's as false (so excludes them).
same_hts_both_fem_inds <- which(same_hts & both_fem)
diff_hts_both_fem_inds <- which(!same_hts & both_fem)
both_male_inds <- which(both_male)
diff_sex_diff_hts_inds <- which(!same_sex & !same_hts)
male_1_fem_2_same_hts_inds <- which(male_1_fem_2 & same_hts)
fem_1_male_2_same_hts_inds <- which(fem_1_male_2 & same_hts)
male_1_fem_2_diff_hts_inds <- which(male_1_fem_2 & !same_hts)
fem_1_male_2_diff_hts_inds <- which(fem_1_male_2 & !same_hts)
male_1_fem_2_inds <- which(male_1_fem_2)
fem_1_male_2_inds <- which(fem_1_male_2)



# Find haplopair probabilities given sexes and an unrelated pair.  Shows up as
# weird frequency table thing when cbind without the as.numeric.
hp_probs_sex_up <- as.numeric(ht_probs_1 * ht_probs_2)



# Find haplopair probabilities given sexes and a parent-offspring pair

# Create vector of NA's for haplopair probabilities given sexes and a
# parent-offspring pair
hp_probs_sex_pop <- rep(NA, n_pairs)

# When both animals are female, and have the same haplotype, the probability is
# the haplotype probability
hp_probs_sex_pop[same_hts_both_fem_inds] <- ht_probs_1[same_hts_both_fem_inds]

# When both animals are female, and have different haplotypes, the probability
# is zero
hp_probs_sex_pop[diff_hts_both_fem_inds] <- 0

# When both animals are male, the probability is the product of the haplotype
# probabilities
hp_probs_sex_pop[both_male_inds] <- ht_probs_1[both_male_inds] *
  ht_probs_2[both_male_inds]

# When one animal is male, and one is female, and their haplotypes are
# different, the probability is half the product of the haplotype probabilities
hp_probs_sex_pop[diff_sex_diff_hts_inds] <- 1/2 *
  ht_probs_1[diff_sex_diff_hts_inds] * ht_probs_2[diff_sex_diff_hts_inds]

# When one animal is male, and one is female, and their haplotypes are the same,
# the probability is the mean of the product of the haplotype probabilities, and
# the probability of the female's haplotype
hp_probs_sex_pop[male_1_fem_2_same_hts_inds] <- 1/2 *
  ht_probs_2[male_1_fem_2_same_hts_inds] *
  (ht_probs_1[male_1_fem_2_same_hts_inds] + 1)
hp_probs_sex_pop[fem_1_male_2_same_hts_inds] <- 1/2 *
  ht_probs_1[fem_1_male_2_same_hts_inds] *
  (ht_probs_2[fem_1_male_2_same_hts_inds] + 1)



# Find haplopair probabilities given sexes and ordered parent-offspring and
# offspring-parent pairs

# Copy from unordered parent-offspring pair probabilities
hp_probs_sex_po <- hp_probs_sex_op <- hp_probs_sex_pop

# When the parent is male the probability is just the product of the haplotype
# probabilities
hp_probs_sex_po[male_1_fem_2_inds] <- 
  ht_probs_1[male_1_fem_2_inds] * ht_probs_2[male_1_fem_2_inds]
hp_probs_sex_op[fem_1_male_2_inds] <- 
  ht_probs_1[fem_1_male_2_inds] * ht_probs_2[fem_1_male_2_inds]

# When the parent is female and the haplotypes are different the probability is
# zero
hp_probs_sex_po[fem_1_male_2_diff_hts_inds] <- 0
hp_probs_sex_op[male_1_fem_2_diff_hts_inds] <- 0

# When the parent is female and the haplotypes are the same the probability is
# the haplotype probability
hp_probs_sex_po[fem_1_male_2_same_hts_inds] <- 
  ht_probs_1[fem_1_male_2_same_hts_inds]
hp_probs_sex_op[male_1_fem_2_same_hts_inds] <- 
  ht_probs_1[male_1_fem_2_same_hts_inds]




# Find haplopair probabilities given sexes and a self-pair

# Find pairs with same sex and haplotype
same_hts_same_sex <- same_hts & same_sex
same_hts_same_sex_inds <- which(same_hts_same_sex)

# Create vector of NA's for haplopair probabilities given sexes and a self-pair
hp_probs_sex_sp <- rep(NA, n_pairs)

# Find haplopair probabilities given sexes and a self-pair.  These are 0 when
# sexes equal but haplotypes not, and NaN when sexes not equal
hp_probs_sex_sp[which(same_sex & !same_hts)] <- 0
hp_probs_sex_sp[which(!same_sex)] <- NaN
hp_probs_sex_sp[same_hts_same_sex_inds] <- ht_probs_1[same_hts_same_sex_inds]



# Find haplopair probabilities given sexes and half-sibs

# Find indices for pairs with equal haplotypes
same_hts_inds <- which(same_hts)

# Find haplopair probabilities given sexes and half-sibs.  Mean of maternal
# case, meaning must be equal, and paternal case, meaning can be anything.  Also
# turns probs tables to numeric for cbinding later.
hp_probs_sex_hsp <- as.numeric(1/2 * ht_probs_1 * ht_probs_2)
hp_probs_sex_hsp[same_hts_inds] <- hp_probs_sex_hsp[same_hts_inds] + 
  1/2 * ht_probs_1[same_hts_inds]



# Find haplopair probabilities given sexes and full-sibs

# Create vector of NA's for haplopair probabilities given sexes and full-sibs
hp_probs_sex_fsp <- rep(NA, n_pairs)

# Find haplopair probabilities given sexes and full-sibs.  Same mum implies must
# be equal
hp_probs_sex_fsp[same_hts_inds] <- ht_probs_1[same_hts_inds]
hp_probs_sex_fsp[which(!same_hts)] <- 0



# Combine log haplopair probabilities given sexes and simple kinships
log_hp_probs_sex_kin <- data.frame(
  unrelated = log(hp_probs_sex_up),
  half_sibs = log(hp_probs_sex_hsp),
  po_op = log(hp_probs_sex_pop),
  po = log(hp_probs_sex_po),
  op = log(hp_probs_sex_op),
  full_sibs = log(hp_probs_sex_fsp),
  self = log(hp_probs_sex_sp)
)

# Remove haplopair probabilities as no longer needed
rm(hp_probs_sex_up, hp_probs_sex_hsp, hp_probs_sex_pop, hp_probs_sex_po,
   hp_probs_sex_op, hp_probs_sex_fsp, hp_probs_sex_sp)



# Find kinship log likelihoods given sexes and haplopairs for simple kinships.
# The sex probabilities are approximated by 0.5, although it is more like 0.55
# for being female.  But these are only relevant to distinguishing self-pairs,
# which is easy with any decent genotype

# Probabilities for missing haplopairs are 1, logs are 0
kin_log_likes_hps_sexes <- log_hp_probs_sex_kin
kin_log_likes_hps_sexes[is.na(log_hp_probs_sex_kin)] <- 0

# Non-self kinships, sex-pairs are independent
sexp_not_na <- !is.na(same_sex)
kin_log_likes_hps_sexes[sexp_not_na, 1:6] <-
  kin_log_likes_hps_sexes[sexp_not_na, 1:6] + log(1/4)

# Self kinships, sexes cannot be different
same_sex_inds <- which(same_sex)
kin_log_likes_hps_sexes[same_sex_inds, 7] <-
  kin_log_likes_hps_sexes[same_sex_inds, 7] + log(1/2)
kin_log_likes_hps_sexes[which(!same_sex), 7] <- -Inf



# Combine kinship log likelihoods given microsatellite genotypes, MtDNA
# haplotypes, and sexes, assuming that the first are conditionally independent
# of the latter two, given kinship. This will give NA if either term is NA, so
# adjust below if not both NA
kin_log_likes_msats_hps_sexes <- kin_log_likes_msats + kin_log_likes_hps_sexes

# Find indices of pairs where exactly one of the terms above is NA
one_na <- is.na(kin_log_likes_msats_hps_sexes)
na_1_not_na_2_inds <- one_na & !is.na(kin_log_likes_hps_sexes)
na_2_not_na_1_inds <- one_na & !is.na(kin_log_likes_msats)

# If only one missing then set to the other
kin_log_likes_msats_hps_sexes[na_1_not_na_2_inds] <- 
  kin_log_likes_hps_sexes[na_1_not_na_2_inds]
kin_log_likes_msats_hps_sexes[na_2_not_na_1_inds] <- 
  kin_log_likes_msats[na_2_not_na_1_inds]



# Save data as R objects
save(list = objects(), file = "rdata/whale_analysis_5.Rdata")