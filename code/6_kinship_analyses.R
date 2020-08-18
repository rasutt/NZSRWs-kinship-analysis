# Find maximum likelihood kinships, numbers of loci with 0, 1, or 2 alleles
# shared, and hsp vs up plods, for each pair of animals

# Find the maximum likelihood kinship for each pair, and the numbers of pairs
# for which each kinship is most likely, given only this data (no
# prior probability)
kinships <- names(kin_log_likes_msats_hps_sexes)
max_like_kin_inds <- apply(kin_log_likes_msats_hps_sexes, 1, which.max)
max_like_kinships <- kinships[max_like_kin_inds]
max_like_kin_tab <- tabulate(max_like_kin_inds, nbins = 7)
names(max_like_kin_tab) <- kinships



# Find numbers of loci where both alleles are shared for all pairs
ns_loci_both_ales_shared <- rowSums(both_ales_shared, na.rm = T)

# Create table of proportions of msat loci with 0, 1, and 2 alleles shared
props_loci_ns_ales_shared <- data.frame(
  data.frame(
    p_0 = ns_loci_both_obs - ns_loci_ales_shared,
    p_1 = (ns_loci_ales_shared - ns_loci_both_ales_shared),
    p_2 = ns_loci_both_ales_shared
  ) / ns_loci_both_obs,
  n_loci = ns_loci_both_obs
)
head(props_loci_ns_ales_shared)



# Find pairs with complete data
pairs_complete <- ns_loci_both_obs == 13 & sexp_not_na & hp_not_na
sum(pairs_complete)



# Find half-siblings versus unrelated pair plods 
hsp_up_plods <- 
  (kin_log_likes_msats_hps_sexes$SOK - kin_log_likes_msats_hps_sexes$UP) / 
  ns_loci_both_obs

# Find PLODs for pairs with complete data
hsp_up_plods_comp <- hsp_up_plods[pairs_complete]

# Find PLODs for possible POPs
poss_pop_hsp_plods <- hsp_up_plods[possible_pop_inds]

# Find PLODs for possible POPs with complete data
poss_pop_hsp_plods_comp <- poss_pop_hsp_plods[
  pairs_complete[possible_pop_inds]
  ]
length(poss_pop_hsp_plods_comp)



# Find half-siblings versus unrelated pair plods 
fsp_up_plods <- 
  (kin_log_likes_msats_hps_sexes$SOKx2 - kin_log_likes_msats_hps_sexes$UP) / 
  ns_loci_both_obs

# Find PLODs for pairs with complete data
fsp_up_plods_comp <- fsp_up_plods[pairs_complete]

# Find PLODs for possible POPs
poss_pop_fsp_plods <- fsp_up_plods[possible_pop_inds]

# Find PLODs for possible POPs with complete data
poss_pop_fsp_plods_comp <- poss_pop_fsp_plods[
  pairs_complete[possible_pop_inds]
  ]
length(poss_pop_fsp_plods_comp)



# Check pairs with largest plods.  Found data problems this way before 
threshold <- 0.5

large_hsp_up_plods_inds <- which(hsp_up_plods > threshold)
large_hsp_up_plods_comp_inds <- which(hsp_up_plods_comp > threshold)

length(large_hsp_up_plods_inds)
length(large_hsp_up_plods_comp_inds)

large_hsp_up_plods <- hsp_up_plods[large_hsp_up_plods_inds]
large_hsp_up_plods_comp <- 
  hsp_up_plods_comp[large_hsp_up_plods_comp_inds]



# Check pairs with smallest plods
threshold <- -0.75

shared_ales[hsp_up_plods < threshold, 1:10]
small_hsp_up_plods <- hsp_up_plods[hsp_up_plods < threshold]
plot(
  density(small_hsp_up_plods)
)
sort(table(ani_pairs_inds[, hsp_up_plods < threshold]))
sum(hsp_up_plods < threshold)
whale_data[389, ]
whale_data[311, ]
?subset

# Next steps:
#   - rewrite this code nicely, 
#   - update the interactive pop connections plot, and 
#   - start inferring possible parent-offspring directionality?

# You can subset matrices and arrays with matrices with columns for coordinates!
# Definitely do this for plods, it'll be so much quicker.  Doesn't matter till
# you're running simulations again though.  Can make a 3d array with the
# matrices in the top left corner and do all loci at once.  Tensors haha!  Can
# probably also do the allele distribution matrix you were thinking about.



# Save data as R objects
save(list = objects(), file = "rdata/whale_analysis_6.Rdata")