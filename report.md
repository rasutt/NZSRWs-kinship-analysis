  - [Microsatellite genotypes, missing data, and
    homozygosity](#microsatellite-genotypes-missing-data-and-homozygosity)
  - [Covariates and haplotypes](#covariates-and-haplotypes)
  - [Shared microsatellite alleles, kinship likelihoods, and
    parent-offspring
    pairs](#shared-microsatellite-alleles-kinship-likelihoods-and-parent-offspring-pairs)

``` r
# Load results of analysis
load(file = "rdata/whale_analysis_6.Rdata")
```

### Microsatellite genotypes, missing data, and homozygosity

``` r
# Find and show proportions of missing data by locus from largest to smallest
barplot(
  colMeans(missing_ales_2),
  las = 2,
  main = "Missing data by locus",
  ylab = "Proportion"
)
```

![](report_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Find the excess observed homozygosity at each locus
props_homo_exp <- sapply(ale_dists, function(dist) sum(dist^2))
props_homo_obs <- colMeans(ales_mat_1 == ales_mat_2, na.rm = T)

# Show observed versus expected homozygosity by locus from largest to smallest
# difference
barplot(
  props_homo_obs - props_homo_exp,
  las = 2,
  main = "Excess observed homozygosity",
  ylab = "Proportion of data"
)
```

![](report_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Locus TR3G1 had the most missing data, and is the only one with
significant excess homozygosity. From Emma’s thesis (pg.84):

“Of the 13 loci retained in the dataset, 12 did not deviate
significantly from Hardy-Weinberg equilibrium and showed no signs of
stutter, allelic dropout or null alleles (Table 2.5). The exception was
TR3G1; this locus had evidence of allelic dropout but was retained as it
was highly informative and allelic dropout was accounted for by
re-amplifying mismatching loci for suspected replicate samples.”

It looks like about half of the apparent homozygotes at this locus may
be due to allelic dropout, about 15% of the data. If we had borderline
cases we could check how this affected them.

### Covariates and haplotypes

The year, location, and field data have been selected randomly when
multiple captures, as Emma just chose the ones with the strongest
genetic results, least missing markers or something, so I’ve removed
them from the dataset for now. I can add them back from multiple
recaptures later if necessary.

Counts for sexes and haplotypes.

``` r
# Show sex and haplotype counts
table(sex_field_hap$SEX)
```

    ## 
    ##   M   F 
    ## 371 449

``` r
table(sex_field_hap$DLP)
```

    ## 
    ##    BakHapA   BakHapB'   BakHapB+    BakHapC    BakHapD    BakHapE    CarHapJ 
    ##        275         96        241         65         91          6          1 
    ##   PatHap04 PatMalHapB   PorHap17       SWPJ 
    ##         20          5          1          1

Samples for which I have not found the original entries with sex or
haplotype data yet.

``` r
# Sort names for animals for which an original entry has not yet been found
as.character(
  whale_data[!(whale_data$Sample.Name %in% sex_field_hap$Sample.Name), 1]
)
```

    ##  [1] "Eau98AI080" "Eau98AI084" "Eau98AI102" "Eau06NZ06"  "Eau08AI017"
    ##  [6] "Eau09NZ06"  "U18-055"    "U13-064"    "U13-065"    "Eau09AI026"
    ## [11] "U13-063"    "Eau09AI159" "Eau98AI150" "U13-073"    "Eau09AI067"
    ## [16] "Eau09AI007" "Eau09AI010" "Eau09AI012" "U18-048"    "U13-066"   
    ## [21] "U13-072"    "Eau05NZ03"

``` r
# Check these with Emma.  Still lots of 2009 animals unfound, where could they be? 
```

There is not much actually missing sex or haplotype data.

``` r
# Check missing sex data
mean(is.na(whale_data$SEX))
```

    ## [1] 0.03276699

``` r
# Check missing haplotype data.
mean(is.na(whale_data$DLP))
```

    ## [1] 0.05339806

### Shared microsatellite alleles, kinship likelihoods, and parent-offspring pairs

Parent-offspring pairs (POPs) share an allele by descent at each locus.
This is not enough to distinguish them from these microsatellite
genotypes, but we can find the pairs that could possibly be POPs
according to them (bearing in mind genotyping error), and check that the
more sophisticated analyses below give the same results.

``` r
# Plot the proportions of loci where pairs are observed to share at least one
# allele
hist(
  ns_loci_ales_shared[ns_loci_both_obs == 13],
  breaks = 100,
  main = "Numbers of loci with at least one shared allele in 
  pairs with complete microsatellite genotypes",
  xlab = "Number of loci",
  ylab = "Number of pairs"
)
```

![](report_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# Check same possible POPs from shared microsatellite alleles and kinship likelihoods
all(which(is.finite(kin_log_likes_msats$PO_OP)) == possible_pop_inds)
```

    ## [1] TRUE

It is possible to have genotypes that share alleles at a locus but which
are still more likely to be unrelated than to be a POP. For example, if
a pair of heterozygous genotypes only shares one allele, the conditional
probability of the second given that they are unrelated is twice the
product of the allele probabilities, and given that they are a POP it is
half the unshared allele probability. So if the probability of the
shared allele is greater than one quarter, the unrelated pair prob is
higher. For SNPs, allele probabilities are often greater than one
quarter, but you can’t have two heterozygous genotypes that only share
one allele.