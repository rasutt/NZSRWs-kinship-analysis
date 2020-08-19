# Distinguish/integrate year, location, field, sex, and haplotype data. The
# first three are currently possibly randomly chosen from among multiple
# recaptures

# Read data from Males All, Females All, and Sheet1 sheets from big multi-sheet
# file, and CI_QC_Unique, and UNIQUE MNZ may 2011 nc sheets from Mainland and
# Campbell Island file.
males_all <- read.csv(
  "spreadsheet_csvs/SRW_database_V2_July2012_NoCalves_repaired - Males All.csv"
)
females_all <- read.csv(
  "spreadsheet_csvs/SRW_database_V2_July2012_NoCalves_repaired - Females All.csv"
)
sheet_1 <- read.csv(
  "spreadsheet_csvs/SRW_database_V2_July2012_NoCalves_repaired - Sheet1.csv"
)
campbell <- read.csv("spreadsheet_csvs/MNZ&CI_For_RAS - CI_QC_Unique.csv")
mainland <- read.csv(
  "spreadsheet_csvs/MNZ&CI_For_RAS - UNIQUE MNZ may 2011 nc.csv"
)

# Standardise column-names over all sheets
names(sheet_1)[5] <- "field.data"
names(sheet_2009)[3:4] <- c("DLP", "SEX")
names(campbell)[c(1:2, 3)] <- names(mainland)[c(1, 4, 3)] <-
  c("Sample.Name", "SEX", "DLP")

# Add empty field data column for Campbell Island animals
campbell$field.data <- ""

# Tidy up sample names for Mainland animals
split_names <- strsplit(mainland$Sample.Name, " |\\(")
mainland$Sample.Name <- sapply(split_names, function(name_vec) name_vec[1])

# Combine sex, field data, and haplotypes from all sheets.  rbind matches by
# name not position
sex_field_hap <- rbind(
  males_all[1:314, c(1, 5:7)], 
  females_all[1:388, c(1, 5:7)], 
  sheet_1[1:57, c(1, 3:5)],
  sheet_2009[1:182, c(1, 3:5)],
  campbell[1:21, c(1:3, 48)],
  mainland[, c(1, 3:4, 6)]
)

# Remove duplicates
sex_field_hap <- sex_field_hap[!duplicated(sex_field_hap$Sample.Name), ]



# Tidy up haplotype, sex, and field data

# Standardise terminology
sex_field_hap$DLP[
  sex_field_hap$DLP %in% c(
    "PorHap4.1", 
    "PorHap4", 
    "PORHAP4",
    "PATHAP4.1",
    "4.1",
    "PORHAP4.1",
    "PATHAP04"
  ) 
  ] <- "PatHap04"
sex_field_hap$DLP[sex_field_hap$DLP %in% c("A", "BAKHAPA")] <- "BakHapA"
sex_field_hap$DLP[sex_field_hap$DLP %in% c("B'", "BAKHAPB'")] <- "BakHapB'"
sex_field_hap$DLP[
  sex_field_hap$DLP %in% c("B+", "B++", "BAKHAPB+", "BAKHAPB++") 
  ] <- "BakHapB+"
sex_field_hap$DLP[sex_field_hap$DLP %in% c("C", "BAKHAPC")] <- "BakHapC"
sex_field_hap$DLP[sex_field_hap$DLP %in% c("D", "BAKHAPD")] <- "BakHapD"
sex_field_hap$DLP[sex_field_hap$DLP %in% c("E", "BAKHAPE")] <- "BakHapE"
sex_field_hap$DLP[sex_field_hap$DLP == "PATMALHAPB"] <- "PatMalHapB"
sex_field_hap$DLP[sex_field_hap$DLP == "PORHAP17"] <- "PorHap17"

sex_field_hap$field.data[sex_field_hap$field.data %in% c("CALF", "calf")] <-
  "Calf"
sex_field_hap$field.data[
  sex_field_hap$field.data %in% c("COW", "cow ", "cow")
  ] <- "Cow"
sex_field_hap$field.data[
  sex_field_hap$field.data %in% c("AD", "Ad", "adult")
  ] <- "Adult"
sex_field_hap$field.data[
  sex_field_hap$field.data %in% c("JUV", "juvenile")
  ] <- "Juvenile"

# Make sex and haplotype columns into factors, excluding missing data
sex_field_hap$SEX <- factor(
  sex_field_hap$SEX,
  exclude = c(NA, "", "N/A", "Faint F")
)
sex_field_hap$DLP <- factor(
  sex_field_hap$DLP,
  exclude = c(NA, "n/a", "N/A", "", "FAIL", "?B+")
)



# Add the haplotype, sex, and field data to the table of all unique animals

# Left-join on sample name.  The rows in the first table without entries in the
# second table get added back at the bottom.
whale_data <- merge(
  whale_data, 
  sex_field_hap, 
  by = "Sample.Name", 
  all.x = T,
  sort = F
)



# Unpack capture locations from sample names of animals

# Find capture location sections of sample names
capture_locations <- substr(whale_data$Sample.Name, 6, 7)

# Convert to factor
capture_locations <- factor(capture_locations, c("AI", "CI", "NZ"))

# Add NA level
capture_locations <- addNA(capture_locations)

# Change labels
levels(capture_locations) <- 
  c("AI", "CI", "MNZ", "Calf") 

# Add to data
whale_data$location <- capture_locations



# Unpack capture years from sample names of animals

# Find capture year sections of sample names
capture_years <- substr(whale_data$Sample.Name, 4, 5)

# Adjust incorrect sections
wrong_cap_years <- capture_years == "-0"
capture_years[wrong_cap_years] <-
  substr(whale_data$Sample.Name[wrong_cap_years], 2, 3)

# Convert to full year numbers
capture_years <- as.integer(capture_years)
capture_years_1900s <- capture_years > 90
capture_years[capture_years_1900s] <- capture_years[capture_years_1900s] + 1900
capture_years[!capture_years_1900s] <- 
  capture_years[!capture_years_1900s] + 2000

# Add to data
whale_data$year <- capture_years

# Reorder columns
whale_data <- whale_data[c(1, 32:31, 29:30, 28, 2:27)]

# Convert updated data to matrices for better access
ales_mat_1 <- as.matrix(whale_data[, 6 + ale_inds_1])
ales_mat_2 <- as.matrix(whale_data[, 6 + ale_inds_2])

# Find missing alleles.  There is one case where the second allele is OOB so
# use missing second-alleles to exclude it with the others
missing_ales_2 <- is.na(ales_mat_2)



# Write data to a csv file
write.csv(whale_data, "combined_whale_data.csv", row.names = F)

# Save data as R objects
save(list = objects(), file = "rdata/whale_analysis_2.Rdata")