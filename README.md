Code to analyse the kinships among a sample of NZ southern right whales using microsatellite genotypes, sex, and MtDNA haplotypes.  

It was written on an old machine with very limited memory and compute, when my main machine broke down without a recent backup, which is why it is the second attempt at the same problem.  It is thus much more carefully developed than the first attempt, including removing large intermediate variables from memory when they are no longer needed.

It's organised into 6 numbered scripts located in the code folder, which each load the previous results and save their own in corresponding files in the rdata folder. 

source_scripts.R in the main folder is used to conveniently call those that are required, in order, when changes are made.

It takes .csv files, derived from excell spreadsheets, as inputs, and returns the likelihoods given to various simple kinships, for each pair of animals, by that data.  

It also does some simple analysis of the possibility that each pair is a parent-offspring pair (POP), based on whether they share an allele at each marker in common.  

It also combines some covariate data including sample year and location.

It saved possible POPs and displayed them in a shiny app, though the code to save that data may have changed since then.

It also writes a single new .csv file with the integrated data for each animal, though that integration is not yet complete. 

I stopped work on it before ISEC 2020 and then Emma gave me the Global SRW RAD genotype data and I started working on that.

Things to do:
- Add to github, maybe combine the report and the readme
- Fix the sex probabilities, they're not 0.5 (mean(whale_data$SEX == "M", na.rm = T) -> 0.455)
- Probably way faster to find possible sex and haplotype combinations and reference them rather than computing them individually, but probably not worth rewriting now.  Quite tricky dealing with missing data.  Really not a big deal, already got sucked into some work that wasn't really worth it.  Maybe just make a note.
- Get the final results sorted
- Would be interesting to compare fsp vs up plods, see if the curves are more informative.  Takes mixed kinships into account more maybe.
- Change likely POPs to possible POPs, depending on final PLOD plots.  
- Update the report and web app.  Don't need too much detail, just the main results/explanation.
- Change all density plots to histograms
- Update folder on gdrive
- Make copy of web app for Global dataset.