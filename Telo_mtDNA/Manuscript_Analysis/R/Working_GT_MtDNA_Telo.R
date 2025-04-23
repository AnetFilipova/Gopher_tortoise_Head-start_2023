#### Exploring the effects of first-year cold-dormancy on telomere length and mitochondrial DNA density in blood cells of head-started Gopher tortoises ####

## Load the necessary libraries ##
library(dplyr) ## Will be used for specifying logical operators
library(stringr) ## Will be used to handle character vectors when modifying text data
library(tidyr) ## Will be used for reshaping the data for easier visualization later on
library(readr) ## Will be used when reading in .csv data
library(lme4) ## Will be used for mixed effect models with random and fixed effects 
library(emmeans) ## Will be used for post-hoc comparisons
library(ggplot2) ## Will be used for making graphs/data visualization 


# Clear memory
rm(list=ls(all = TRUE))


#Process EDIT all the runs for EACH Plate of samples sepearately
#Use the ...Quantification Cq Results.csv files
#Copy and paste the block of code and repeat on each plate of samples.
#You can use find replace for the dataset name, Plate1 for Plate2


## Plate 1##
# Import .csv files for each run for a particular plate of samples 
Plate1_MPX <- read.csv("GT_Plate1_Multiplex_12_11_2024_Quantification_Cq_Results.csv")
dim (Plate1_MPX)

Plate1_Telo <- read.csv("GT_Plate1_Telomeres_12_11_2024_Quantification_Cq_Results.csv")
dim (Plate1_Telo)

# Concatenate data across runs for the same plate of samples
Plate1<- rbind(Plate1_MPX, Plate1_Telo)
dim(Plate1)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate1$PlateID <- "Plate1"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate1$Target[Plate1$Fluor == "VIC"] <- "scnag"
Plate1$Target[Plate1$Fluor == "FAM"] <- "mtdna"
Plate1$Target[Plate1$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate1$Diff_AVG_Cq <- abs(Plate1$Cq - Plate1$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.4 threshold
Plate1$Flag_outlier <- ifelse(Plate1$Target %in% c("mtdna", "scnag") & Plate1$Diff_AVG_Cq > 0.4, "yes", "no")

# Identify outliers for telomeres based on >0.4 threshold
Plate1$Flag_outlier <- ifelse(Plate1$Target == "telomeres" & Plate1$Diff_AVG_Cq > 0.4, "yes", Plate1$Flag_outlier)

# Report and examine high Cq samples in the first round
HighCq <- Plate1[Plate1$Flag_outlier == "yes", c("Well", "Sample", "Fluor", "Diff_AVG_Cq")]
print(paste("Number of rows with high Cq in the first round:", nrow(HighCq)))
print(HighCq)

# Remove outliers identified in the first round
Plate1 <- Plate1[Plate1$Flag_outlier != "yes", ]

# Verify the dimensions after removing outliers in the first round
print(dim(Plate1))

############## Look for outliers - 2nd round ##############

# Recalculate Cq Mean and Starting Quantity Mean based on remaining data
Plate1 <- Plate1 %>%
  group_by(Sample, Target) %>%
  mutate(
    Cq.Mean2 = mean(Cq),
    Sq.Mean2 = mean(Starting.Quantity..SQ.)
  ) %>%
  ungroup()

# Calculate the new difference from the updated mean Cq
Plate1$Diff_AVG_Cq_2 <- abs(Plate1$Cq - Plate1$Cq.Mean2)

# Identify outliers in the second round based on >0.4 threshold for telomeres
Plate1$Flag_outlier_2 <- ifelse(Plate1$Target == "telomeres" & Plate1$Diff_AVG_Cq_2 > 0.4, "yes", "no")

# Report and examine high Cq samples in the second round for telomeres
HighCq_2 <- Plate1[Plate1$Target == "telomeres" & Plate1$Diff_AVG_Cq_2 > 0.4, c("Well", "Sample", "Fluor", "Diff_AVG_Cq_2")]
print(paste("Number of rows with high Cq in the second round for telomeres:", nrow(HighCq_2)))
print(HighCq_2)

# Remove outliers identified in the second round for telomeres
Plate1 <- Plate1[!(Plate1$Target == "telomeres" & Plate1$Flag_outlier_2 == "yes"), ]

# Verify the dimensions after removing outliers in the second round
print(dim(Plate1))

### Remove samples that do not have at least two rows
Plate1 <- Plate1 %>%
  group_by(Sample) %>%
  filter(n() >= 2) %>%
  ungroup()

# Final dimensions after all filtering steps
print(dim(Plate1))

######### Remove negative controls and standards
# rows that have "NEG", "POS" in column "Sample" and remove rows with "STD" in Sample "Content"
Plate1 <- Plate1 %>%
  filter(!str_detect(Sample, "STD")) %>%
  filter(!str_detect(Sample, "NEG")) %>%
  filter(!str_detect(Sample, "POS")) 
dim(Plate1)

######### Subset dataset based on the value in "Target" column
unique_targets <- unique(Plate1$Target)

# Create a list to store the subset dataframes
subset_dfs <- list()

# Loop through each unique value in 'Target', subset the dataframe, and store in subset_dfs
for (target_value in unique_targets) {
  subset_df <- subset(Plate1, Target == target_value)
  subset_dfs[[target_value]] <- subset_df
}


####### Now subset_dfs is a list where each element is a dataframe containing rows for each unique 'Target' value
Plate1_SCNAG<-print(subset_dfs[["scnag"]])
Plate1_SCNAG<-Plate1_SCNAG[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean",  "Flag_outlier", "SQ.Mean")]
Plate1_SCNAG <- Plate1_SCNAG %>% 
  rename(Target_SCNAG = Target, Cq_SCNAG = Cq, Cq.Mean_SCNAG = Cq.Mean, Flag_outlier_SCNAG=Flag_outlier, SQ.Mean_SCNAG=SQ.Mean)

Plate1_mtDNA<-print(subset_dfs[["mtdna"]])
Plate1_mtDNA<-Plate1_mtDNA[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate1_mtDNA <- Plate1_mtDNA %>% 
  rename(Target_mtDNA = Target, Cq_mtDNA = Cq, Cq.Mean_mtDNA = Cq.Mean, Flag_outlier_mtDNA = Flag_outlier, SQ.Mean_mtDNA = SQ.Mean)

Plate1_Telomeres<-print(subset_dfs[["telomeres"]])
Plate1_Telomeres<-Plate1_Telomeres[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate1_Telomeres <- Plate1_Telomeres %>%
  rename(Target_Telomeres = Target, Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)
  #rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)
  

#############  Make a final MPX dataset for by merging the Target datasets horizontally, in rows. 
Plate1_FinalMPX <- merge(Plate1_SCNAG, Plate1_mtDNA, by = c("PlateID", "Well", "Sample"))

############# Normalize mtDNA
# Add a column called mtDNA, and calculate the normalized value 
Plate1_FinalMPX$mtDNA <- (Plate1_FinalMPX$SQ.Mean_mtDNA / Plate1_FinalMPX$SQ.Mean_SCNAG)

# Recalculate mean across the replicates
Plate1_FinalMPX <- Plate1_FinalMPX %>%
  group_by(Sample) %>%
  mutate(
    mtDNA.Mean = mean(mtDNA)) %>%
  ungroup()

############## Merge final MPX with Telomeres
## Reduce datasets to single row per individual containing only the columns we want.
Plate1_FinalMPX <- distinct(Plate1_FinalMPX, PlateID, Sample, SQ.Mean_SCNAG, Cq.Mean_SCNAG, SQ.Mean_mtDNA, mtDNA.Mean)
Plate1_FinalTelo <- distinct(Plate1_Telomeres, PlateID, Sample, SQ.Mean_Telomeres, Cq.Mean_Telomeres)

## Merge the files horizontally
Plate1_FinalData <- merge(Plate1_FinalMPX, Plate1_FinalTelo, by = c("PlateID", "Sample"))

# Normalize Telomeres
Plate1_FinalData <- Plate1_FinalData %>% mutate(Telomeres.per.cell = SQ.Mean_Telomeres / SQ.Mean_SCNAG)

# Assuming Plate1_FinalData has been prepared as per your previous steps

# Aggregate to get one row per sample, taking mean of normalized mtDNA and telomeres
Plate1_FinalData <- Plate1_FinalData %>%
  group_by(Sample) %>%
  summarize(
    SQ.Mean_SCNAG = mean(SQ.Mean_SCNAG),
    SQ.Mean_mtDNA = mean(SQ.Mean_mtDNA),
    mtDNA.Mean = mean(mtDNA.Mean),
    Cq.Mean_SCNAG = mean(Cq.Mean_SCNAG),
    SQ.Mean_Telomeres = mean(SQ.Mean_Telomeres),
    Cq.Mean_Telomeres = mean(Cq.Mean_Telomeres),
    Telomeres.per.cell = mean(Telomeres.per.cell)
  ) %>%
  ungroup()

# Write the final data file for this plate
write.csv(file = "Plate1_FinalData.csv", Plate1_FinalData, row.names = FALSE)

#########################################################

#### Plate 2 ####

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate2_MPX <- read.csv("GT_Plate2_Multiplex_12_12_2024_Quantification_Cq_Results.csv")
dim (Plate2_MPX)

Plate2_Telo <- read.csv("GT_Plate2_Telomeres_12_12_2024_Quantification_Cq_Results.csv")
dim (Plate2_Telo)

# Concatenate data across runs for the same plate of samples
Plate2<- rbind(Plate2_MPX, Plate2_Telo)
dim(Plate2)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate2$PlateID <- "Plate2"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate2$Target[Plate2$Fluor == "VIC"] <- "scnag"
Plate2$Target[Plate2$Fluor == "FAM"] <- "mtdna"
Plate2$Target[Plate2$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate2$Diff_AVG_Cq <- abs(Plate2$Cq - Plate2$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.4 threshold
Plate2$Flag_outlier <- ifelse(Plate2$Target %in% c("mtdna", "scnag") & Plate2$Diff_AVG_Cq > 0.4, "yes", "no")

# Identify outliers for telomeres based on >0.4 threshold
Plate2$Flag_outlier <- ifelse(Plate2$Target == "telomeres" & Plate2$Diff_AVG_Cq > 0.4, "yes", Plate2$Flag_outlier)

# Report and examine high Cq samples in the first round
HighCq <- Plate2[Plate2$Flag_outlier == "yes", c("Well", "Sample", "Fluor", "Diff_AVG_Cq")]
print(paste("Number of rows with high Cq in the first round:", nrow(HighCq)))
print(HighCq)

# Remove outliers identified in the first round
Plate2 <- Plate2[Plate2$Flag_outlier != "yes", ]

# Verify the dimensions after removing outliers in the first round
print(dim(Plate2))

############## Look for outliers - 2nd round ##############

# Recalculate Cq Mean and Starting Quantity Mean based on remaining data
Plate2 <- Plate2 %>%
  group_by(Sample, Target) %>%
  mutate(
    Cq.Mean2 = mean(Cq),
    Sq.Mean2 = mean(Starting.Quantity..SQ.)
  ) %>%
  ungroup()

# Calculate the new difference from the updated mean Cq
Plate2$Diff_AVG_Cq_2 <- abs(Plate2$Cq - Plate2$Cq.Mean2)

# Identify outliers in the second round based on >0.4 threshold for telomeres
Plate2$Flag_outlier_2 <- ifelse(Plate2$Target == "telomeres" & Plate2$Diff_AVG_Cq_2 > 0.4, "yes", "no")

# Report and examine high Cq samples in the second round for telomeres
HighCq_2 <- Plate2[Plate2$Target == "telomeres" & Plate2$Diff_AVG_Cq_2 > 0.4, c("Well", "Sample", "Fluor", "Diff_AVG_Cq_2")]
print(paste("Number of rows with high Cq in the second round for telomeres:", nrow(HighCq_2)))
print(HighCq_2)

# Remove outliers identified in the second round for telomeres
Plate2 <- Plate2[!(Plate2$Target == "telomeres" & Plate2$Flag_outlier_2 == "yes"), ]

# Verify the dimensions after removing outliers in the second round
print(dim(Plate2))

### Remove samples that do not have at least two rows
Plate2 <- Plate2 %>%
  group_by(Sample) %>%
  filter(n() >= 2) %>%
  ungroup()

# Final dimensions after all filtering steps
print(dim(Plate2))

######### Remove negative controls and standards
# rows that have "NEG", "POS" in column "Sample" and remove rows with "STD" in Sample "Content"
Plate2 <- Plate2 %>%
  filter(!str_detect(Sample, "STD")) %>%
  filter(!str_detect(Sample, "NEG")) 
dim(Plate2)

######### Subset dataset based on the value in "Target" column
unique_targets <- unique(Plate2$Target)

# Create a list to store the subset dataframes
subset_dfs <- list()

# Loop through each unique value in 'Target', subset the dataframe, and store in subset_dfs

for (target_value in unique_targets) {
  subset_df <- subset(Plate2, Target == target_value)
  subset_dfs[[target_value]] <- subset_df
}

#############
#############Now subset_dfs is a list where each element is a dataframe containing rows for each unique 'Target' value
Plate2_SCNAG<-print(subset_dfs[["scnag"]])
Plate2_SCNAG<-Plate2_SCNAG[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean",  "Flag_outlier", "SQ.Mean")]
Plate2_SCNAG <- Plate2_SCNAG %>% 
  rename(Target_SCNAG = Target, Cq_SCNAG = Cq, Cq.Mean_SCNAG = Cq.Mean, Flag_outlier_SCNAG=Flag_outlier, SQ.Mean_SCNAG=SQ.Mean)

Plate2_mtDNA<-print(subset_dfs[["mtdna"]])
Plate2_mtDNA<-Plate2_mtDNA[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate2_mtDNA <- Plate2_mtDNA %>% 
  rename(Target_mtDNA = Target, Cq_mtDNA = Cq, Cq.Mean_mtDNA = Cq.Mean, Flag_outlier_mtDNA = Flag_outlier, SQ.Mean_mtDNA = SQ.Mean)

Plate2_Telomeres<-print(subset_dfs[["telomeres"]])
Plate2_Telomeres<-Plate2_Telomeres[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate2_Telomeres <- Plate2_Telomeres %>% 
  rename(Target_Telomeres = Target, Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)
  #rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)

###############. Make a final MPX dataset for by merging the Target datasets horizontally, in rows. 
Plate2_FinalMPX <- merge(Plate2_SCNAG, Plate2_mtDNA, by = c("PlateID", "Well", "Sample"))

############# Normalize mtDNA
# Add a column called mtDNA, and calculate the normalized value 
Plate2_FinalMPX$mtDNA <- (Plate2_FinalMPX$SQ.Mean_mtDNA / Plate2_FinalMPX$SQ.Mean_SCNAG)

# Recalculate mean across the replicates
Plate2_FinalMPX <- Plate2_FinalMPX %>%
  group_by(Sample) %>%
  mutate(
    mtDNA.Mean = mean(mtDNA)) %>%
  ungroup()

############## Merge final MPX with Telomeres
## Reduce datasets to single row per individual containing only the columns we want.
Plate2_FinalMPX <- distinct(Plate2_FinalMPX, PlateID, Sample, SQ.Mean_SCNAG, Cq.Mean_SCNAG, SQ.Mean_mtDNA, mtDNA.Mean)
Plate2_FinalTelo <- distinct(Plate2_Telomeres, PlateID, Sample, SQ.Mean_Telomeres, Cq.Mean_Telomeres)

## Merge the files horizontally
Plate2_FinalData <- merge(Plate2_FinalMPX, Plate2_FinalTelo, by = c("PlateID", "Sample"))

# Normalize Telomeres
Plate2_FinalData <- Plate2_FinalData %>% mutate(Telomeres.per.cell = SQ.Mean_Telomeres / SQ.Mean_SCNAG)

# Assuming Plate2_FinalData has been prepared as per your previous steps

# Aggregate to get one row per sample, taking mean of normalized mtDNA and telomeres
Plate2_FinalData <- Plate2_FinalData %>%
  group_by(Sample) %>%
  summarize(
    SQ.Mean_SCNAG = mean(SQ.Mean_SCNAG),
    SQ.Mean_mtDNA = mean(SQ.Mean_mtDNA),
    mtDNA.Mean = mean(mtDNA.Mean),
    Cq.Mean_SCNAG = mean(Cq.Mean_SCNAG),
    SQ.Mean_Telomeres = mean(SQ.Mean_Telomeres),
    Cq.Mean_Telomeres = mean(Cq.Mean_Telomeres),
    Telomeres.per.cell = mean(Telomeres.per.cell)
  ) %>%
  ungroup()

# Write the final data file for this plate
write.csv(file = "Plate2_FinalData.csv", Plate2_FinalData, row.names = FALSE)

#########################################################

## ## Plate 3##

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate3_MPX <- read.csv("GT_Plate3_Multiplex_12_16_2024_Quantification_Cq_Results.csv")
dim (Plate3_MPX)

Plate3_Telo <- read.csv("GT_Plate3_Telomeres_12_16_2024_Quantification_Cq_Results.csv")
dim (Plate3_Telo)

# Concatenate data across runs for the same plate of samples
Plate3<- rbind(Plate3_MPX, Plate3_Telo)
dim(Plate3)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate3$PlateID <- "Plate3"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate3$Target[Plate3$Fluor == "VIC"] <- "scnag"
Plate3$Target[Plate3$Fluor == "FAM"] <- "mtdna"
Plate3$Target[Plate3$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate3$Diff_AVG_Cq <- abs(Plate3$Cq - Plate3$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.4 threshold
Plate3$Flag_outlier <- ifelse(Plate3$Target %in% c("mtdna", "scnag") & Plate3$Diff_AVG_Cq > 0.4, "yes", "no")

# Identify outliers for telomeres based on >0.4 threshold
Plate3$Flag_outlier <- ifelse(Plate3$Target == "telomeres" & Plate3$Diff_AVG_Cq > 0.4, "yes", Plate3$Flag_outlier)

# Report and examine high Cq samples in the first round
HighCq <- Plate3[Plate3$Flag_outlier == "yes", c("Well", "Sample", "Fluor", "Diff_AVG_Cq")]
print(paste("Number of rows with high Cq in the first round:", nrow(HighCq)))
print(HighCq)

# Remove outliers identified in the first round
Plate3 <- Plate3[Plate3$Flag_outlier != "yes", ]

# Verify the dimensions after removing outliers in the first round
print(dim(Plate3))

############## Look for outliers - 2nd round ##############

# Recalculate Cq Mean and Starting Quantity Mean based on remaining data
Plate3 <- Plate3 %>%
  group_by(Sample, Target) %>%
  mutate(
    Cq.Mean2 = mean(Cq),
    Sq.Mean2 = mean(Starting.Quantity..SQ.)
  ) %>%
  ungroup()

# Calculate the new difference from the updated mean Cq
Plate3$Diff_AVG_Cq_2 <- abs(Plate3$Cq - Plate3$Cq.Mean2)

# Identify outliers in the second round based on >0.4 threshold for telomeres
Plate3$Flag_outlier_2 <- ifelse(Plate3$Target == "telomeres" & Plate3$Diff_AVG_Cq_2 > 0.4, "yes", "no")

# Report and examine high Cq samples in the second round for telomeres
HighCq_2 <- Plate3[Plate3$Target == "telomeres" & Plate3$Diff_AVG_Cq_2 > 0.4, c("Well", "Sample", "Fluor", "Diff_AVG_Cq_2")]
print(paste("Number of rows with high Cq in the second round for telomeres:", nrow(HighCq_2)))
print(HighCq_2)

# Remove outliers identified in the second round for telomeres
Plate3 <- Plate3[!(Plate3$Target == "telomeres" & Plate3$Flag_outlier_2 == "yes"), ]

# Verify the dimensions after removing outliers in the second round
print(dim(Plate3))

### Remove samples that do not have at least two rows
Plate3 <- Plate3 %>%
  group_by(Sample) %>%
  filter(n() >= 2) %>%
  ungroup()

# Final dimensions after all filtering steps
print(dim(Plate3))

######### Remove negative controls and standards
# rows that have "NEG", "POS" in column "Sample" and remove rows with "STD" in Sample "Content"
Plate3 <- Plate3 %>%
  filter(!str_detect(Sample, "STD")) %>%
  filter(!str_detect(Sample, "NEG")) 
dim(Plate3)

######### Subset dataset based on the value in "Target" column
unique_targets <- unique(Plate3$Target)
# Create a list to store the subset dataframes
subset_dfs <- list()
# Loop through each unique value in 'Target', subset the dataframe, and store in subset_dfs
for (target_value in unique_targets) {
  subset_df <- subset(Plate3, Target == target_value)
  subset_dfs[[target_value]] <- subset_df
}

############# Now subset_dfs is a list where each element is a dataframe containing rows for each unique 'Target' value
Plate3_SCNAG<-print(subset_dfs[["scnag"]])
Plate3_SCNAG<-Plate3_SCNAG[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean",  "Flag_outlier", "SQ.Mean")]
Plate3_SCNAG <- Plate3_SCNAG %>% 
  rename(Target_SCNAG = Target, Cq_SCNAG = Cq, Cq.Mean_SCNAG = Cq.Mean, Flag_outlier_SCNAG=Flag_outlier, SQ.Mean_SCNAG=SQ.Mean)

Plate3_mtDNA<-print(subset_dfs[["mtdna"]])
Plate3_mtDNA<-Plate3_mtDNA[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate3_mtDNA <- Plate3_mtDNA %>% 
  rename(Target_mtDNA = Target, Cq_mtDNA = Cq, Cq.Mean_mtDNA = Cq.Mean, Flag_outlier_mtDNA = Flag_outlier, SQ.Mean_mtDNA = SQ.Mean)

Plate3_Telomeres<-print(subset_dfs[["telomeres"]])
Plate3_Telomeres<-Plate3_Telomeres[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate3_Telomeres <- Plate3_Telomeres %>% 
  rename(Target_Telomeres = Target, Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)
  #rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)

###############. Make a final MPX dataset for by merging the Target datasets horizontally, in rows. 
Plate3_FinalMPX <- merge(Plate3_SCNAG, Plate3_mtDNA, by = c("PlateID", "Well", "Sample"))

############# Normalize mtDNA
# Add a column called mtDNA, and calculate the normalized value 
Plate3_FinalMPX$mtDNA <- (Plate3_FinalMPX$SQ.Mean_mtDNA / Plate3_FinalMPX$SQ.Mean_SCNAG)

# Recalculate mean across the replicates
Plate3_FinalMPX <- Plate3_FinalMPX %>%
  group_by(Sample) %>%
  mutate(
    mtDNA.Mean = mean(mtDNA)) %>%
  ungroup()

############## Merge final MPX with Telomeres
## Reduce datasets to single row per individual containing only the columns we want.
Plate3_FinalMPX <- distinct(Plate3_FinalMPX, PlateID, Sample, SQ.Mean_SCNAG, Cq.Mean_SCNAG, SQ.Mean_mtDNA, mtDNA.Mean)
Plate3_FinalTelo <- distinct(Plate3_Telomeres, PlateID, Sample, SQ.Mean_Telomeres, Cq.Mean_Telomeres)

## Merge the files horizontally
Plate3_FinalData <- merge(Plate3_FinalMPX, Plate3_FinalTelo, by = c("PlateID", "Sample"))

# Normalize Telomeres
Plate3_FinalData <- Plate3_FinalData %>% mutate(Telomeres.per.cell = SQ.Mean_Telomeres / SQ.Mean_SCNAG)

# Assuming Plate3_FinalData has been prepared as per your previous steps

# Aggregate to get one row per sample, taking mean of normalized mtDNA and telomeres
Plate3_FinalData <- Plate3_FinalData %>%
  group_by(Sample) %>%
  summarize(
    SQ.Mean_SCNAG = mean(SQ.Mean_SCNAG),
    SQ.Mean_mtDNA = mean(SQ.Mean_mtDNA),
    mtDNA.Mean = mean(mtDNA.Mean),
    Cq.Mean_SCNAG = mean(Cq.Mean_SCNAG),
    SQ.Mean_Telomeres = mean(SQ.Mean_Telomeres),
    Cq.Mean_Telomeres = mean(Cq.Mean_Telomeres),
    Telomeres.per.cell = mean(Telomeres.per.cell)
  ) %>%
  ungroup()

# Write the final data file for this plate
write.csv(file = "Plate3_FinalData.csv", Plate3_FinalData, row.names = FALSE)

#########################################################

## ## Plate 4##

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate4_MPX <- read.csv("GT_Plate4_Multiplex_12_17_2024_Quantification_Cq_Results.csv")
dim (Plate4_MPX)

Plate4_Telo <- read.csv("GT_Plate4_Telomeres_12_17_2024_Quantification_Cq_Results.csv")
dim (Plate4_Telo)

# Concatenate data across runs for the same plate of samples
Plate4<- rbind(Plate4_MPX, Plate4_Telo)
dim(Plate4)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate4$PlateID <- "Plate4"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate4$Target[Plate4$Fluor == "VIC"] <- "scnag"
Plate4$Target[Plate4$Fluor == "FAM"] <- "mtdna"
Plate4$Target[Plate4$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate4$Diff_AVG_Cq <- abs(Plate4$Cq - Plate4$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.4 threshold
Plate4$Flag_outlier <- ifelse(Plate4$Target %in% c("mtdna", "scnag") & Plate4$Diff_AVG_Cq > 0.4, "yes", "no")

# Identify outliers for telomeres based on >0.4 threshold
Plate4$Flag_outlier <- ifelse(Plate4$Target == "telomeres" & Plate4$Diff_AVG_Cq > 0.4, "yes", Plate4$Flag_outlier)

# Report and examine high Cq samples in the first round
HighCq <- Plate4[Plate4$Flag_outlier == "yes", c("Well", "Sample", "Fluor", "Diff_AVG_Cq")]
print(paste("Number of rows with high Cq in the first round:", nrow(HighCq)))
print(HighCq)

# Remove outliers identified in the first round
Plate4 <- Plate4[Plate4$Flag_outlier != "yes", ]

# Verify the dimensions after removing outliers in the first round
print(dim(Plate4))

############## Look for outliers - 2nd round ##############

# Recalculate Cq Mean and Starting Quantity Mean based on remaining data
Plate4 <- Plate4 %>%
  group_by(Sample, Target) %>%
  mutate(
    Cq.Mean2 = mean(Cq),
    Sq.Mean2 = mean(Starting.Quantity..SQ.)
  ) %>%
  ungroup()

# Calculate the new difference from the updated mean Cq
Plate4$Diff_AVG_Cq_2 <- abs(Plate4$Cq - Plate4$Cq.Mean2)

# Identify outliers in the second round based on >0.4 threshold for telomeres
Plate4$Flag_outlier_2 <- ifelse(Plate4$Target == "telomeres" & Plate4$Diff_AVG_Cq_2 > 0.4, "yes", "no")

# Report and examine high Cq samples in the second round for telomeres
HighCq_2 <- Plate4[Plate4$Target == "telomeres" & Plate4$Diff_AVG_Cq_2 > 0.4, c("Well", "Sample", "Fluor", "Diff_AVG_Cq_2")]
print(paste("Number of rows with high Cq in the second round for telomeres:", nrow(HighCq_2)))
print(HighCq_2)

# Remove outliers identified in the second round for telomeres
Plate4 <- Plate4[!(Plate4$Target == "telomeres" & Plate4$Flag_outlier_2 == "yes"), ]

# Verify the dimensions after removing outliers in the second round
print(dim(Plate4))

### Remove samples that do not have at least two rows
Plate4 <- Plate4 %>%
  group_by(Sample) %>%
  filter(n() >= 2) %>%
  ungroup()

# Final dimensions after all filtering steps
print(dim(Plate4))

######### Remove negative controls and standards
# rows that have "NEG", "POS" in column "Sample" and remove rows with "STD" in Sample "Content"
Plate4 <- Plate4 %>%
  filter(!str_detect(Sample, "STD")) %>%
  filter(!str_detect(Sample, "NEG")) 
dim(Plate4)

######### Subset dataset based on the value in "Target" column
unique_targets <- unique(Plate4$Target)
# Create a list to store the subset dataframes
subset_dfs <- list()

# Loop through each unique value in 'Target', subset the dataframe, and store in subset_dfs
for (target_value in unique_targets) {
  subset_df <- subset(Plate4, Target == target_value)
  subset_dfs[[target_value]] <- subset_df
}

############# Now subset_dfs is a list where each element is a dataframe containing rows for each unique 'Target' value
Plate4_SCNAG<-print(subset_dfs[["scnag"]])
Plate4_SCNAG<-Plate4_SCNAG[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean",  "Flag_outlier", "SQ.Mean")]
Plate4_SCNAG <- Plate4_SCNAG %>% 
  rename(Target_SCNAG = Target, Cq_SCNAG = Cq, Cq.Mean_SCNAG = Cq.Mean, Flag_outlier_SCNAG=Flag_outlier, SQ.Mean_SCNAG=SQ.Mean)

Plate4_mtDNA<-print(subset_dfs[["mtdna"]])
Plate4_mtDNA<-Plate4_mtDNA[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate4_mtDNA <- Plate4_mtDNA %>% 
  rename(Target_mtDNA = Target, Cq_mtDNA = Cq, Cq.Mean_mtDNA = Cq.Mean, Flag_outlier_mtDNA = Flag_outlier, SQ.Mean_mtDNA = SQ.Mean)

Plate4_Telomeres<-print(subset_dfs[["telomeres"]])
Plate4_Telomeres<-Plate4_Telomeres[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate4_Telomeres <- Plate4_Telomeres %>%
  rename(Target_Telomeres = Target, Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)
  #rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)

############### Make a final MPX dataset for by merging the Target datasets horizontally, in rows. 
Plate4_FinalMPX <- merge(Plate4_SCNAG, Plate4_mtDNA, by = c("PlateID", "Well", "Sample"))

############# Normalize mtDNA
# Add a column called mtDNA, and calculate the normalized value 
Plate4_FinalMPX$mtDNA <- (Plate4_FinalMPX$SQ.Mean_mtDNA / Plate4_FinalMPX$SQ.Mean_SCNAG)

# Recalculate mean across the replicates
Plate4_FinalMPX <- Plate4_FinalMPX %>%
  group_by(Sample) %>%
  mutate(
    mtDNA.Mean = mean(mtDNA)) %>%
  ungroup()

############## Merge final MPX with Telomeres
## Reduce datasets to single row per individual containing only the columns we want.
Plate4_FinalMPX <- distinct(Plate4_FinalMPX, PlateID, Sample, SQ.Mean_SCNAG, Cq.Mean_SCNAG, SQ.Mean_mtDNA, mtDNA.Mean)
Plate4_FinalTelo <- distinct(Plate4_Telomeres, PlateID, Sample, SQ.Mean_Telomeres, Cq.Mean_Telomeres)

## Merge the files horizontally
Plate4_FinalData <- merge(Plate4_FinalMPX, Plate4_FinalTelo, by = c("PlateID", "Sample"))

# Normalize Telomeres
Plate4_FinalData <- Plate4_FinalData %>% mutate(Telomeres.per.cell = SQ.Mean_Telomeres / SQ.Mean_SCNAG)

# Assuming Plate4_FinalData has been prepared as per your previous steps

# Aggregate to get one row per sample, taking mean of normalized mtDNA and telomeres
Plate4_FinalData <- Plate4_FinalData %>%
  group_by(Sample) %>%
  summarize(
    SQ.Mean_SCNAG = mean(SQ.Mean_SCNAG),
    SQ.Mean_mtDNA = mean(SQ.Mean_mtDNA),
    mtDNA.Mean = mean(mtDNA.Mean),
    Cq.Mean_SCNAG = mean(Cq.Mean_SCNAG),
    SQ.Mean_Telomeres = mean(SQ.Mean_Telomeres),
    Cq.Mean_Telomeres = mean(Cq.Mean_Telomeres),
    Telomeres.per.cell = mean(Telomeres.per.cell)
  ) %>%
  ungroup()

# Write the final data file for this plate
write.csv(file = "Plate4_FinalData.csv", Plate4_FinalData, row.names = FALSE)

############################################
## ## Plate 5##

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate5_MPX <- read.csv("GT_Plate5_Multiplex_12_17_2024_Quantification_Cq_Results.csv")
dim (Plate5_MPX)

Plate5_Telo <- read.csv("GT_Plate5_Telomeres_12_17_2024_Quantification_Cq_Results.csv")
dim (Plate5_Telo)

# Concatenate data across runs for the same plate of samples
Plate5<- rbind(Plate5_MPX, Plate5_Telo)
dim(Plate5)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate5$PlateID <- "Plate5"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate5$Target[Plate5$Fluor == "VIC"] <- "scnag"
Plate5$Target[Plate5$Fluor == "FAM"] <- "mtdna"
Plate5$Target[Plate5$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate5$Diff_AVG_Cq <- abs(Plate5$Cq - Plate5$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.4 threshold
Plate5$Flag_outlier <- ifelse(Plate5$Target %in% c("mtdna", "scnag") & Plate5$Diff_AVG_Cq > 0.4, "yes", "no")

# Identify outliers for telomeres based on >0.4 threshold
Plate5$Flag_outlier <- ifelse(Plate5$Target == "telomeres" & Plate5$Diff_AVG_Cq > 0.4, "yes", Plate5$Flag_outlier)

# Report and examine high Cq samples in the first round
HighCq <- Plate5[Plate5$Flag_outlier == "yes", c("Well", "Sample", "Fluor", "Diff_AVG_Cq")]
print(paste("Number of rows with high Cq in the first round:", nrow(HighCq)))
print(HighCq)

# Remove outliers identified in the first round
Plate5 <- Plate5[Plate5$Flag_outlier != "yes", ]

# Verify the dimensions after removing outliers in the first round
print(dim(Plate5))

############## Look for outliers - 2nd round ##############

# Recalculate Cq Mean and Starting Quantity Mean based on remaining data
Plate5 <- Plate5 %>%
  group_by(Sample, Target) %>%
  mutate(
    Cq.Mean2 = mean(Cq),
    Sq.Mean2 = mean(Starting.Quantity..SQ.)
  ) %>%
  ungroup()

# Calculate the new difference from the updated mean Cq
Plate5$Diff_AVG_Cq_2 <- abs(Plate5$Cq - Plate5$Cq.Mean2)

# Identify outliers in the second round based on >0.4 threshold for telomeres
Plate5$Flag_outlier_2 <- ifelse(Plate5$Target == "telomeres" & Plate5$Diff_AVG_Cq_2 > 0.4, "yes", "no")

# Report and examine high Cq samples in the second round for telomeres
HighCq_2 <- Plate5[Plate5$Target == "telomeres" & Plate5$Diff_AVG_Cq_2 > 0.4, c("Well", "Sample", "Fluor", "Diff_AVG_Cq_2")]
print(paste("Number of rows with high Cq in the second round for telomeres:", nrow(HighCq_2)))
print(HighCq_2)

# Remove outliers identified in the second round for telomeres
Plate5 <- Plate5[!(Plate5$Target == "telomeres" & Plate5$Flag_outlier_2 == "yes"), ]

# Verify the dimensions after removing outliers in the second round
print(dim(Plate5))

### Remove samples that do not have at least two rows
Plate5 <- Plate5 %>%
  group_by(Sample) %>%
  filter(n() >= 2) %>%
  ungroup()

# Final dimensions after all filtering steps
print(dim(Plate5))

######### Remove negative controls and standards
# rows that have "NEG", "POS" in column "Sample" and remove rows with "STD" in Sample "Content"
Plate5 <- Plate5 %>%
  filter(!str_detect(Sample, "STD")) %>%
  filter(!str_detect(Sample, "NEG"))
dim(Plate5)

######### Subset dataset based on the value in "Target" column
unique_targets <- unique(Plate5$Target)

# Create a list to store the subset dataframes
subset_dfs <- list()

# Loop through each unique value in 'Target', subset the dataframe, and store in subset_dfs
for (target_value in unique_targets) {
  subset_df <- subset(Plate5, Target == target_value)
  subset_dfs[[target_value]] <- subset_df
}

############# Now subset_dfs is a list where each element is a dataframe containing rows for each unique 'Target' value
Plate5_SCNAG<-print(subset_dfs[["scnag"]])
Plate5_SCNAG<-Plate5_SCNAG[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean",  "Flag_outlier", "SQ.Mean")]
Plate5_SCNAG <- Plate5_SCNAG %>% 
  rename(Target_SCNAG = Target, Cq_SCNAG = Cq, Cq.Mean_SCNAG = Cq.Mean, Flag_outlier_SCNAG=Flag_outlier, SQ.Mean_SCNAG=SQ.Mean)

Plate5_mtDNA<-print(subset_dfs[["mtdna"]])
Plate5_mtDNA<-Plate5_mtDNA[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate5_mtDNA <- Plate5_mtDNA %>% 
  rename(Target_mtDNA = Target, Cq_mtDNA = Cq, Cq.Mean_mtDNA = Cq.Mean, Flag_outlier_mtDNA = Flag_outlier, SQ.Mean_mtDNA = SQ.Mean)

Plate5_Telomeres<-print(subset_dfs[["telomeres"]])
Plate5_Telomeres<-Plate5_Telomeres[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate5_Telomeres <- Plate5_Telomeres %>% 
  rename(Target_Telomeres = Target, Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)
  #rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)

###############. Make a final MPX dataset for by merging the Target datasets horizontally, in rows. 
Plate5_FinalMPX <- merge(Plate5_SCNAG, Plate5_mtDNA, by = c("PlateID", "Well", "Sample"))

############# Normalize mtDNA
# Add a column called mtDNA, and calculate the normalized value 
Plate5_FinalMPX$mtDNA <- (Plate5_FinalMPX$SQ.Mean_mtDNA / Plate5_FinalMPX$SQ.Mean_SCNAG)

# Recalculate mean across the replicates
Plate5_FinalMPX <- Plate5_FinalMPX %>%
  group_by(Sample) %>%
  mutate(
    mtDNA.Mean = mean(mtDNA)) %>%
  ungroup()

############## Merge final MPX with Telomeres
## Reduce datasets to single row per individual containing only the columns we want.
Plate5_FinalMPX <- distinct(Plate5_FinalMPX, PlateID, Sample, SQ.Mean_SCNAG, Cq.Mean_SCNAG, SQ.Mean_mtDNA, mtDNA.Mean)
Plate5_FinalTelo <- distinct(Plate5_Telomeres, PlateID, Sample, SQ.Mean_Telomeres, Cq.Mean_Telomeres)

## Merge the files horizontally
Plate5_FinalData <- merge(Plate5_FinalMPX, Plate5_FinalTelo, by = c("PlateID", "Sample"))

# Normalize Telomeres
Plate5_FinalData <- Plate5_FinalData %>% mutate(Telomeres.per.cell = SQ.Mean_Telomeres / SQ.Mean_SCNAG)

# Assuming Plate5_FinalData has been prepared as per your previous steps

# Aggregate to get one row per sample, taking mean of normalized mtDNA and telomeres
Plate5_FinalData <- Plate5_FinalData %>%
  group_by(Sample) %>%
  summarize(
    SQ.Mean_SCNAG = mean(SQ.Mean_SCNAG),
    SQ.Mean_mtDNA = mean(SQ.Mean_mtDNA),
    mtDNA.Mean = mean(mtDNA.Mean),
    Cq.Mean_SCNAG = mean(Cq.Mean_SCNAG),
    SQ.Mean_Telomeres = mean(SQ.Mean_Telomeres),
    Cq.Mean_Telomeres = mean(Cq.Mean_Telomeres),
    Telomeres.per.cell = mean(Telomeres.per.cell)
  ) %>%
  ungroup()

# Write the final data file for this plate
write.csv(file = "Plate5_FinalData.csv", Plate5_FinalData, row.names = FALSE)

#######################################
## ## Plate 6##

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate6_MPX <- read.csv("GT_Plate6_Multiplex_12_18_2024_Quantification_Cq_Results.csv")
dim (Plate6_MPX)

Plate6_MPX <- Plate6_MPX[, -c(1)]

Plate6_Telo <- read.csv("GT_Plate6_Telomeres_12_18_2024_Quantification_Cq_Results.csv")
dim (Plate6_Telo)

# Concatenate data across runs for the same plate of samples
Plate6<- rbind(Plate6_MPX, Plate6_Telo)
dim(Plate6)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate6$PlateID <- "Plate6"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate6$Target[Plate6$Fluor == "VIC"] <- "scnag"
Plate6$Target[Plate6$Fluor == "FAM"] <- "mtdna"
Plate6$Target[Plate6$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate6$Diff_AVG_Cq <- abs(Plate6$Cq - Plate6$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.4 threshold
Plate6$Flag_outlier <- ifelse(Plate6$Target %in% c("mtdna", "scnag") & Plate6$Diff_AVG_Cq > 0.4, "yes", "no")

# Identify outliers for telomeres based on >0.4 threshold
Plate6$Flag_outlier <- ifelse(Plate6$Target == "telomeres" & Plate6$Diff_AVG_Cq > 0.4, "yes", Plate6$Flag_outlier)

# Report and examine high Cq samples in the first round
HighCq <- Plate6[Plate6$Flag_outlier == "yes", c("Well", "Sample", "Fluor", "Diff_AVG_Cq")]
print(paste("Number of rows with high Cq in the first round:", nrow(HighCq)))
print(HighCq)

# Remove outliers identified in the first round
Plate6 <- Plate6[Plate6$Flag_outlier != "yes", ]

# Verify the dimensions after removing outliers in the first round
print(dim(Plate6))

############## Look for outliers - 2nd round ##############

# Recalculate Cq Mean and Starting Quantity Mean based on remaining data
Plate6 <- Plate6 %>%
  group_by(Sample, Target) %>%
  mutate(
    Cq.Mean2 = mean(Cq),
    Sq.Mean2 = mean(Starting.Quantity..SQ.)
  ) %>%
  ungroup()

# Calculate the new difference from the updated mean Cq
Plate6$Diff_AVG_Cq_2 <- abs(Plate6$Cq - Plate6$Cq.Mean2)

# Identify outliers in the second round based on >0.4 threshold for telomeres
Plate6$Flag_outlier_2 <- ifelse(Plate6$Target == "telomeres" & Plate6$Diff_AVG_Cq_2 > 0.4, "yes", "no")

# Report and examine high Cq samples in the second round for telomeres
HighCq_2 <- Plate6[Plate6$Target == "telomeres" & Plate6$Diff_AVG_Cq_2 > 0.4, c("Well", "Sample", "Fluor", "Diff_AVG_Cq_2")]
print(paste("Number of rows with high Cq in the second round for telomeres:", nrow(HighCq_2)))
print(HighCq_2)

# Remove outliers identified in the second round for telomeres
Plate6 <- Plate6[!(Plate6$Target == "telomeres" & Plate6$Flag_outlier_2 == "yes"), ]

# Verify the dimensions after removing outliers in the second round
print(dim(Plate6))

### Remove samples that do not have at least two rows
Plate6 <- Plate6 %>%
  group_by(Sample) %>%
  filter(n() >= 2) %>%
  ungroup()

# Final dimensions after all filtering steps
print(dim(Plate6))

######### Remove negative controls and standards
# rows that have "NEG", "POS" in column "Sample" and remove rows with "STD" in Sample "Content"
Plate6 <- Plate6 %>%
  filter(!str_detect(Sample, "STD")) %>%
  filter(!str_detect(Sample, "NEG")) 
dim(Plate6)

######### Subset dataset based on the value in "Target" column
unique_targets <- unique(Plate6$Target)

# Create a list to store the subset dataframes
subset_dfs <- list()

# Loop through each unique value in 'Target', subset the dataframe, and store in subset_dfs
for (target_value in unique_targets) {
  subset_df <- subset(Plate6, Target == target_value)
  subset_dfs[[target_value]] <- subset_df
}

############# Now subset_dfs is a list where each element is a dataframe containing rows for each unique 'Target' value
Plate6_SCNAG<-print(subset_dfs[["scnag"]])
Plate6_SCNAG<-Plate6_SCNAG[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean",  "Flag_outlier", "SQ.Mean")]
Plate6_SCNAG <- Plate6_SCNAG %>% 
  rename(Target_SCNAG = Target, Cq_SCNAG = Cq, Cq.Mean_SCNAG = Cq.Mean, Flag_outlier_SCNAG=Flag_outlier, SQ.Mean_SCNAG=SQ.Mean)

Plate6_mtDNA<-print(subset_dfs[["mtdna"]])
Plate6_mtDNA<-Plate6_mtDNA[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate6_mtDNA <- Plate6_mtDNA %>% 
  rename(Target_mtDNA = Target, Cq_mtDNA = Cq, Cq.Mean_mtDNA = Cq.Mean, Flag_outlier_mtDNA = Flag_outlier, SQ.Mean_mtDNA = SQ.Mean)

Plate6_Telomeres<-print(subset_dfs[["telomeres"]])
Plate6_Telomeres<-Plate6_Telomeres[ ,c("PlateID", "Well", "Sample", "Target", "Cq", "Cq.Mean", "Flag_outlier", "SQ.Mean")]
Plate6_Telomeres <- Plate6_Telomeres %>% 
  rename(Target_Telomeres = Target, Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)
  #rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)

############### Make a final MPX dataset for by merging the Target datasets horizontally, in rows. 
Plate6_FinalMPX <- merge(Plate6_SCNAG, Plate6_mtDNA, by = c("PlateID", "Well", "Sample"))

############# Normalize mtDNA
# Add a column called mtDNA, and calculate the normalized value 
Plate6_FinalMPX$mtDNA <- (Plate6_FinalMPX$SQ.Mean_mtDNA / Plate6_FinalMPX$SQ.Mean_SCNAG)

# Recalculate mean across the replicates
Plate6_FinalMPX <- Plate6_FinalMPX %>%
  group_by(Sample) %>%
  mutate(
    mtDNA.Mean = mean(mtDNA)) %>%
  ungroup()

############## Merge final MPX with Telomeres
## Reduce datasets to single row per individual containing only the columns we want.
Plate6_FinalMPX <- distinct(Plate6_FinalMPX, PlateID, Sample, SQ.Mean_SCNAG, Cq.Mean_SCNAG, SQ.Mean_mtDNA, mtDNA.Mean)
Plate6_FinalTelo <- distinct(Plate6_Telomeres, PlateID, Sample, SQ.Mean_Telomeres, Cq.Mean_Telomeres)

## Merge the files horizontally
Plate6_FinalData <- merge(Plate6_FinalMPX, Plate6_FinalTelo, by = c("PlateID", "Sample"))

# Normalize Telomeres
Plate6_FinalData <- Plate6_FinalData %>% mutate(Telomeres.per.cell = SQ.Mean_Telomeres / SQ.Mean_SCNAG)

# Assuming Plate6_FinalData has been prepared as per your previous steps

# Aggregate to get one row per sample, taking mean of normalized mtDNA and telomeres
Plate6_FinalData <- Plate6_FinalData %>%
  group_by(Sample) %>%
  summarize(
    SQ.Mean_SCNAG = mean(SQ.Mean_SCNAG),
    SQ.Mean_mtDNA = mean(SQ.Mean_mtDNA),
    mtDNA.Mean = mean(mtDNA.Mean),
    Cq.Mean_SCNAG = mean(Cq.Mean_SCNAG),
    SQ.Mean_Telomeres = mean(SQ.Mean_Telomeres),
    Cq.Mean_Telomeres = mean(Cq.Mean_Telomeres),
    Telomeres.per.cell = mean(Telomeres.per.cell)
  ) %>%
  ungroup()

# Write the final data file for this plate
write.csv(file = "Plate6_FinalData.csv", Plate6_FinalData, row.names = FALSE)

# Merge all FinalData together, for each plate, by concatenating the rows

# Load each plate's final data
Plate1_FinalData <- read.csv("Plate1_FinalData.csv")
Plate2_FinalData <- read.csv("Plate2_FinalData.csv")
Plate3_FinalData <- read.csv("Plate3_FinalData.csv")
Plate4_FinalData <- read.csv("Plate4_FinalData.csv")
Plate5_FinalData <- read.csv("Plate5_FinalData.csv")
Plate6_FinalData <- read.csv("Plate6_FinalData.csv")

# Concatenate the datasets
qPCR_FinalData <- rbind(Plate1_FinalData, Plate2_FinalData, Plate3_FinalData, Plate4_FinalData, Plate5_FinalData, Plate6_FinalData)

# Verify the dimensions of the combined data
print(dim(qPCR_FinalData))

# Merge Final Data with Trait MetaData for your individuals
# Import .csv files for each run for a particular plate of samples 
Trait <- read.csv("Trait_MetaData.csv")
dim (Trait)
dim(qPCR_FinalData)

# Merge both datasets
FinalData <- merge(qPCR_FinalData, Trait, by = c("Sample"))

# Save the merged data to a final CSV file
write.csv(FinalData, "GT_FinalData.csv", row.names = FALSE)

# Optional: Print the dimensions of the final merged data
print(dim(FinalData))

#################### Calculating the percent change in mtDNA and Telomere length ####################

######### Prepare the data ##########

# Load the data
data <- read_csv("GT_FinalData.csv")

# Creating a new column called Tortoise_ID by extracting the individual ID from the column Sample; sub() is a base R function that replaces the first match of a pattern in a string.
# "_.*" is a regular expression that means: “an underscore (_) followed by anything (.*)”
# "" means: replace that part with nothing (i.e., remove it).

data <- data %>%
  mutate(Tortoise_ID = sub("_.*", "", Sample))

# Arrange by Tortoise_ID and Time_Point (A to D)
data <- data %>%
  arrange(Tortoise_ID, Time_Point)

# Filter data for Interval A to B
interval_A_to_B <- data %>%
  filter(Time_Point %in% c("A", "B"))

# Filter data for Interval B to C
interval_B_to_C <- data %>%
  filter(Time_Point %in% c("B", "C"))

# Filter data for Interval C to D
interval_C_to_D <- data %>%
  filter(Time_Point %in% c("C", "D"))

# Calculate percent changes from one timepoint to the next, i.e. from timepoint A to B
interval_A_to_B <- interval_A_to_B %>%
  group_by(Tortoise_ID) %>%
  mutate(
    mtDNA_percent_change = (mtDNA.Mean - lag(mtDNA.Mean)) / lag(mtDNA.Mean) * 100,
    Telomeres_percent_change = (Telomeres.per.cell - lag(Telomeres.per.cell)) / lag(Telomeres.per.cell) * 100,
    Transition = paste0(lag(Time_Point), "_to_", Time_Point)
  ) %>%
  filter(!is.na(mtDNA_percent_change), !is.na(Telomeres_percent_change)) %>%
  ungroup()

# Calculate percent changes from one timepoint to the next, i.e. from timepoint B to C
interval_B_to_C <- interval_B_to_C %>%
  group_by(Tortoise_ID) %>%
  mutate(
    mtDNA_percent_change = (mtDNA.Mean - lag(mtDNA.Mean)) / lag(mtDNA.Mean) * 100,
    Telomeres_percent_change = (Telomeres.per.cell - lag(Telomeres.per.cell)) / lag(Telomeres.per.cell) * 100,
    Transition = paste0(lag(Time_Point), "_to_", Time_Point)
  ) %>%
  filter(!is.na(mtDNA_percent_change), !is.na(Telomeres_percent_change)) %>%
  ungroup()

# Calculate percent changes from one timepoint to the next, i.e. from timepoint C to D
interval_C_to_D <- interval_C_to_D %>%
  group_by(Tortoise_ID) %>%
  mutate(
    mtDNA_percent_change = (mtDNA.Mean - lag(mtDNA.Mean)) / lag(mtDNA.Mean) * 100,
    Telomeres_percent_change = (Telomeres.per.cell - lag(Telomeres.per.cell)) / lag(Telomeres.per.cell) * 100,
    Transition = paste0(lag(Time_Point), "_to_", Time_Point)
  ) %>%
  filter(!is.na(mtDNA_percent_change), !is.na(Telomeres_percent_change)) %>%
  ungroup()

# Merge all intervals together
merged_intervals <-rbind(interval_A_to_B, interval_B_to_C, interval_C_to_D)

# Now pivot to wide format
percent_change_wide <- merged_intervals %>%
  select(Tortoise_ID, Transition, mtDNA_percent_change, Telomeres_percent_change) %>% # selecting those columns
  pivot_wider(
    names_from = Transition,
    values_from = c(mtDNA_percent_change, Telomeres_percent_change),
    names_glue = "{.value}_{Transition}" #.value refers to the column names from the values_from argument (in this case: mtDNA_percent_change and Telomeres_percent_change)
    #Transition refers to the unique values in the Transition column, which will become part of the new column names.
  )

# Produce a .csv file with one row per Tortoise_ID and each column represents a change in either mtDNA or telomere length from one timepoint to the next
write_csv(percent_change_wide, "GT_PercentChange.csv")


###### Joining the Percent Change file with the Growth Rate file to run the mixed-effect model ######

# Removing Not_viable individuals or those with deviations from normal growth
growth_data <- read_csv("GT_GrowthData.csv") %>%
  filter(!(Tortoise_ID %in% c("Not_Viable", "GT2023_N05.03", "GT2023_N06.01", "GT2023_N05.06", "GT2023_N15.04"))) %>%
  drop_na()  # removes any remaining rows with NA

# Split the ID into two parts, pad both to two digits, and paste back together
percent_change_wide <- read_csv("GT_PercentChange.csv") %>%
  mutate(
    Tortoise_ID = sprintf("GT2023_N%02d.%02d",
                          as.integer(sub("\\..*", "", Tortoise_ID)),  # part before the decimal
                          as.integer(sub(".*\\.", "", Tortoise_ID))   # part after the decimal
    )
  )

# Now merge the two datasets by Tortoise_ID
merged_data <- left_join(growth_data, percent_change_wide, by = "Tortoise_ID")

# Produce a merged .csv file containing growth rate and percent change in mtDNA and telomere length from one timepoint to the next
write_csv(merged_data, "GT_Growth_Telo_mtDNA.csv")

######################### DATA ANALYSIS ############################

################# Telomere length ##########################

# Running three linear regression models at each interval with Percent change in Telomere length as a response variable and growth rate and treatment as independent variables

### Hypothesis: Hatchling tortoises experiencing early life fast growth at constant warm temperature will have shorter telomere length in blood cells 3 months post-dormancy compared to animals that experienced cold dormancy?

merged_data <- read.csv("GT_Growth_Telo_mtDNA.csv", na.strings = "NA")
str(merged_data)

########## Correlation for telomeres Winter (between Before and During) ########
# Subset to remove NAs for interval A to B
subset_telo_ab <- subset(merged_data, !is.na(Telomeres_percent_change_A_to_B))

# Overall correlation
cor.test(subset_telo_ab$Telomeres_percent_change_A_to_B, subset_telo_ab$Growth_rate_During, method = "spearman")

# Correlation by treatment
by(subset_telo_ab, subset_telo_ab$Treatment, function(df) {
  cor.test(df$Telomeres_percent_change_A_to_B, df$Growth_rate_During, method = "spearman")
})

########### Correlation for telomeres Early Spring (between During and 3 Weeks Post)
# Subset to remove NAs for interval B to C
subset_telo_bc <- subset(merged_data, !is.na(Telomeres_percent_change_B_to_C))

# Overall correlation
cor.test(subset_telo_bc$Telomeres_percent_change_B_to_C, subset_telo_bc$Growth_rate_3_Weeks_Post, method = "spearman")

# Correlation by treatment
by(subset_telo_bc, subset_telo_bc$Treatment, function(df) {
  cor.test(df$Telomeres_percent_change_B_to_C, df$Growth_rate_3_Weeks_Post, method = "spearman")
})

######### Correlation for telomeres Late Spring (Between 3 Weeks Post and 3 Months Post)
# Subset to remove NAs for interval C to D
subset_telo_cd <- subset(merged_data, !is.na(Telomeres_percent_change_C_to_D))

# Overall correlation
cor.test(subset_telo_cd$Telomeres_percent_change_C_to_D, subset_telo_cd$Growth_rate_3_Months_Post, method = "spearman")

# Correlation by treatment
by(subset_telo_cd, subset_telo_cd$Treatment, function(df) {
  cor.test(df$Telomeres_percent_change_C_to_D, df$Growth_rate_3_Months_Post, method = "spearman")
})

# Running a linear-regression for A to B interval  
model_telo_ab <- lmer(Telomeres_percent_change_A_to_B ~ Growth_rate_During + Treatment + (1|Nest_ID), data = subset_telo_ab)
summary(model_telo_ab)

emmeans(model_telo_ab, pairwise ~ Treatment)

# Running a linear-regression for B to C interval  
model_telo_bc <- lmer(Telomeres_percent_change_B_to_C ~ Growth_rate_3_Weeks_Post + Treatment + (1|Nest_ID), data = subset_telo_bc)
summary(model_telo_bc)

emmeans(model_telo_bc, pairwise ~ Treatment)

# Running a linear-regression for C to D interval  
model_telo_cd <- lmer(Telomeres_percent_change_C_to_D ~ Growth_rate_3_Months_Post + Treatment + (1|Nest_ID), data = subset_telo_cd)
summary(model_telo_cd)

emmeans(model_telo_cd, pairwise ~ Treatment)


################# Mitochondrial DNA density ##########################

# Running three mixed-effect models at each interval with Percent change in mtDNA as a response variable and Treatment as a fixed effect, and Nest_ID as a random effect

### Hypothesis: Cold dormancy will suppress metabolic phenotypes including lower mitochondrial copy number and the suppressive effect will persist 3 weeks and 3 months after dormancy.

merged_data <- read.csv("GT_Growth_Telo_mtDNA.csv", na.strings = "NA")
str(merged_data)

########## Correlation for mtDNA in Winter (between Before and During) ########
# Subset to remove NAs for interval A to B
subset_mtDNA_ab <- subset(merged_data, !is.na(mtDNA_percent_change_A_to_B))

# Overall correlation
cor.test(subset_mtDNA_ab$mtDNA_percent_change_A_to_B, subset_mtDNA_ab$Growth_rate_During, method = "spearman")

# Correlation by treatment
by(subset_mtDNA_ab, subset_mtDNA_ab$Treatment, function(df) {
  cor.test(df$mtDNA_percent_change_A_to_B, df$Growth_rate_During, method = "spearman")
})

########### Correlation for mtDNA Early Spring (between During and 3 Weeks Post)
# Subset to remove NAs for interval B to C
subset_mtDNA_bc <- subset(merged_data, !is.na(mtDNA_percent_change_B_to_C))

# Overall correlation
cor.test(subset_mtDNA_bc$mtDNA_percent_change_B_to_C, subset_mtDNA_bc$Growth_rate_3_Weeks_Post, method = "spearman")

# Correlation by treatment
by(subset_mtDNA_bc, subset_mtDNA_bc$Treatment, function(df) {
  cor.test(df$mtDNA_percent_change_B_to_C, df$Growth_rate_3_Weeks_Post, method = "spearman")
})

######### Correlation for mtDNA Late Spring (Between 3 Weeks Post and 3 Months Post)
# Subset to remove NAs for interval C to D
subset_mtDNA_cd <- subset(merged_data, !is.na(mtDNA_percent_change_C_to_D))

# Overall correlation
cor.test(subset_mtDNA_cd$mtDNA_percent_change_C_to_D, subset_mtDNA_cd$Growth_rate_3_Months_Post, method = "spearman")

# Correlation by treatment
by(subset_mtDNA_cd, subset_mtDNA_cd$Treatment, function(df) {
  cor.test(df$mtDNA_percent_change_C_to_D, df$Growth_rate_3_Months_Post, method = "spearman")
})

# Running a linear-regression for A to B interval  
model_mtDNA_ab <- lmer(mtDNA_percent_change_A_to_B ~ Treatment + (1|Nest_ID), data = subset_mtDNA_ab)
summary(model_mtDNA_ab)

emmeans(model_mtDNA_ab, pairwise ~ Treatment)

# Running a linear-regression for B to C interval  
model_mtDNA_bc <- lmer(mtDNA_percent_change_B_to_C ~ Treatment + (1|Nest_ID), data = subset_mtDNA_bc)
summary(model_mtDNA_bc)

emmeans(model_mtDNA_bc, pairwise ~ Treatment)

# Running a linear-regression for C to D interval  
model_mtDNA_cd <- lmer(mtDNA_percent_change_C_to_D ~ Treatment + (1|Nest_ID), data = subset_mtDNA_cd)
summary(model_mtDNA_cd)

emmeans(model_mtDNA_cd, pairwise ~ Treatment)


##### Plotting the results for mtDNA  #####

# Tidy the data for plotting
plot_df <- bind_rows(
  subset_mtDNA_ab %>%
    select(Treatment, mtDNA_percent_change_A_to_B) %>%
    mutate(Timepoint = "Winter", PercentChange = mtDNA_percent_change_A_to_B),
  
  subset_mtDNA_bc %>%
    select(Treatment, mtDNA_percent_change_B_to_C) %>%
    mutate(Timepoint = "Early Spring", PercentChange = mtDNA_percent_change_B_to_C),
  
  subset_mtDNA_cd %>%
    select(Treatment, mtDNA_percent_change_C_to_D) %>%
    mutate(Timepoint = "Late Spring", PercentChange = mtDNA_percent_change_C_to_D)
) %>% 
  select(Treatment, Timepoint, PercentChange)

# Set timepoint order
plot_df$Timepoint <- factor(plot_df$Timepoint, levels = c("Winter", "Early Spring", "Late Spring"))

# Summarize for plotting
summary_df <- plot_df %>%
  group_by(Treatment, Timepoint) %>%
  summarise(
    mean = mean(PercentChange, na.rm = TRUE),
    se = sd(PercentChange, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


# Define color-blind friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Define timepoint labels
timepoint_labels <- c(
  "Winter" = "Winter",
  "Early Spring" = "Early Spring",
  "Late Spring" = "Late Spring"
)

# Create the plot
mtDNA_plot <- ggplot(plot_df, aes(x = Timepoint, y = PercentChange, color = Treatment)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    position = position_dodge(0.8),
    fill = "white",
    alpha = 0.5,
    linewidth = 0.5
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
    size = 2,
    alpha = 0.4
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2,
    aes(group = Treatment),
    position = position_dodge(width = 0.8)
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = Treatment),
    linewidth = 0.6,
    position = position_dodge(width = 0.8)
  ) +
  scale_color_manual(
    values = cbbPalette[c(6, 7)],
    name = "Treatment",
    labels = c("Cold-Dormancy", "Constant-Warmth")
  ) +
  scale_x_discrete(labels = timepoint_labels) +
  scale_y_continuous(breaks = seq(-100, 200, by = 20)) +  # Adjust range if needed
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 30, b = 15)),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "right",
    legend.title.align = 0.5
  ) +
  labs(
    x = "Timepoint",
    y = "Percent Change in mtDNA density"
  )

mtDNA_plot

# Save the plot as PNG for final figure production
ggsave(mtDNA_plot, file = "mtDNA_density.png", width = 12, height = 6, dpi = 600)


##### Plotting the results for telomeres #####

# Tidy the data for plotting
plot_df <- bind_rows(
  subset_telo_ab %>%
    select(Treatment, Telomeres_percent_change_A_to_B) %>%
    mutate(Timepoint = "Winter", PercentChange = Telomeres_percent_change_A_to_B),
  
  subset_telo_bc %>%
    select(Treatment, Telomeres_percent_change_B_to_C) %>%
    mutate(Timepoint = "Early Spring", PercentChange = Telomeres_percent_change_B_to_C),
  
  subset_telo_cd %>%
    select(Treatment, Telomeres_percent_change_C_to_D) %>%
    mutate(Timepoint = "Late Spring", PercentChange = Telomeres_percent_change_C_to_D)
) %>% 
  select(Treatment, Timepoint, PercentChange)

# Set timepoint order
plot_df$Timepoint <- factor(plot_df$Timepoint, levels = c("Winter", "Early Spring", "Late Spring"))

# Summarize for plotting
summary_df <- plot_df %>%
  group_by(Treatment, Timepoint) %>%
  summarise(
    mean = mean(PercentChange, na.rm = TRUE),
    se = sd(PercentChange, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Define color-blind friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Define timepoint labels
timepoint_labels <- c(
  "Winter" = "Winter",
  "Early Spring" = "Early Spring",
  "Late Spring" = "Late Spring"
)

# Create the plot
telomeres_plot <- ggplot(plot_df, aes(x = Timepoint, y = PercentChange, color = Treatment)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    position = position_dodge(0.8),
    fill = "white",
    alpha = 0.5,
    linewidth = 0.5
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
    size = 2,
    alpha = 0.4
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2,
    aes(group = Treatment),
    position = position_dodge(width = 0.8)
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = Treatment),
    linewidth = 0.6,
    position = position_dodge(width = 0.8)
  ) +
  scale_color_manual(
    values = cbbPalette[c(6, 7)],
    name = "Treatment",
    labels = c("Cold-Dormancy", "Constant-Warmth")
  ) +
  scale_x_discrete(labels = timepoint_labels) +
  scale_y_continuous(breaks = seq(-100, 500, by = 50)) +  # Adjust range if needed
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 30, b = 15)),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "right",
    legend.title.align = 0.5
  ) +
  labs(
    x = "Timepoint",
    y = "Percent Change in telomeres"
  )

telomeres_plot

# Save the plot as PNG for final figure production
ggsave(telomeres_plot, file = "Telomere_length.png", width = 12, height = 6, dpi = 600)
