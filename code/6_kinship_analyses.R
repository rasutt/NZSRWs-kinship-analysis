# Find maximum likelihood kinships, numbers of microsatellite loci with 0, 1, or
# 2 alleles shared, and half/full-siblings versus unrelated-pair PLODs, for each
# pair of animals

# Find the kinships analysed
kinships <- names(kin_log_likes_msats_hps_sexes)

# Find the index of the maximum likelihood kinship for each pair based on each
# set of data
max_like_kin_inds_msats <- apply(kin_log_likes_msats, 1, which.max)
max_like_kin_inds_hps_sexes <- apply(kin_log_likes_hps_sexes, 1, which.max)
max_like_kin_inds <- apply(kin_log_likes_msats_hps_sexes, 1, which.max)

# Find the maximum likelihood kinship for each pair based on each set of data
max_like_kinships_msats <- kinships[max_like_kin_inds_msats]
max_like_kinships_hps_sexes <- kinships[max_like_kin_inds_hps_sexes]
max_like_kinships <- kinships[max_like_kin_inds]

# Find the numbers of pairs for which each kinship is most likely given
# this data
max_like_kin_tab <- tabulate(max_like_kin_inds, nbins = 7)
names(max_like_kin_tab) <- kinships



# Find numbers of loci where both alleles are shared for all pairs
ns_loci_both_ales_shared <- rowSums(both_ales_shared, na.rm = T)

# Create table of proportions of microsatellite loci with 0, 1, and 2 alleles
# shared
props_loci_ns_ales_shared <- round(
  data.frame(
    data.frame(
      p_0 = ns_loci_both_obs - ns_loci_ales_shared,
      p_1 = (ns_loci_ales_shared - ns_loci_both_ales_shared),
      p_2 = ns_loci_both_ales_shared
    ) / ns_loci_both_obs,
    n_loci = ns_loci_both_obs
  ), 
  2
)



# Find pairs with complete data (microsatellites, haplotypes, and sexes)
pairs_complete <- ns_loci_both_obs == 13 & !is.na(same_hts) & sexp_not_na 



# Find half-siblings versus unrelated pair plods 
hsp_up_plods <- 
  (kin_log_likes_msats_hps_sexes$SOK - kin_log_likes_msats_hps_sexes$UP) / 
  ns_loci_both_obs

# Find PLODs for pairs with complete data
hsp_up_plods_comp <- hsp_up_plods[pairs_complete]



# Save data as R objects
save(list = objects(), file = "rdata/whale_analysis_6.Rdata")