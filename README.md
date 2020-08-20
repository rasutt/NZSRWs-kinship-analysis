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

The code could probably be sped up a lot by finding all possible geno-, haplo-, and sex-pair probabilities just once, and then just indexing them for each observed pair.  

It could also probably be simplified by setting NA probabilities to 1 and/or their logs to 0.  It only needs NAs sometimes where adding probabilities over cases, so it can be set to 1 or 0 at the end of such calculations.

These might not be worth fixing unless/until I end up running simulations or something.

Things to do:
- Update the report, readme, and web app depending on results.  Don't need too much detail, just the main results/explanation.  Maybe combine the report and the readme?
- Make copy of web app for Global dataset.