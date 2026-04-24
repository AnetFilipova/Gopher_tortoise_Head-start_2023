#### Exploring the effects of first-year cold-dormancy on the change in telomere length and mitochondrial DNA density in blood cells of head-started Gopher tortoises ####

## Load the necessary libraries ##
library(lme4)        
library(lmerTest)    
library(emmeans)     
library(car)        
library(scales)      
library(ggplot2)     
library(readr)       
library(tidyr)     
library(stringr)     
library(dplyr)
library(patchwork)   
library(cowplot)    

## Use dplyr versions of these functions ##
rename <- dplyr::rename
mutate <- dplyr::mutate

# Clear memory
rm(list=ls(all = TRUE))


#Process EDIT all the runs for EACH Plate of samples sepearately
#Use the ...Quantification Cq Results.csv files
#Copy and paste the block of code and repeat on each plate of samples.
#You can use find replace for the dataset name, Plate1 for Plate2


## Plate 1##
# Import .csv files for each run for a particular plate of samples 
Plate1_MPX <- read.csv("Telo_mtDNA/Data/GT_Plate1_Multiplex_12_11_2024_Quantification_Cq_Results.csv")
dim (Plate1_MPX)

Plate1_Telo <- read.csv("Telo_mtDNA/Data/GT_Plate1_Telomeres_12_11_2024_Quantification_Cq_Results.csv")
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
write.csv(file = "Telo_mtDNA/Data/Plate1_FinalData.csv", Plate1_FinalData, row.names = FALSE)

#########################################################

#### Plate 2 ####

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate2_MPX <- read.csv("Telo_mtDNA/Data/GT_Plate2_Multiplex_12_12_2024_Quantification_Cq_Results.csv")
dim (Plate2_MPX)

Plate2_Telo <- read.csv("Telo_mtDNA/Data/GT_Plate2_Telomeres_12_12_2024_Quantification_Cq_Results.csv")
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
write.csv(file = "Telo_mtDNA/Data/Plate2_FinalData.csv", Plate2_FinalData, row.names = FALSE)

#########################################################

## ## Plate 3##

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate3_MPX <- read.csv("Telo_mtDNA/Data/GT_Plate3_Multiplex_12_16_2024_Quantification_Cq_Results.csv")
dim (Plate3_MPX)

Plate3_Telo <- read.csv("Telo_mtDNA/Data/GT_Plate3_Telomeres_12_16_2024_Quantification_Cq_Results.csv")
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
write.csv(file = "Telo_mtDNA/Data/Plate3_FinalData.csv", Plate3_FinalData, row.names = FALSE)

#########################################################

## ## Plate 4##

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate4_MPX <- read.csv("Telo_mtDNA/Data/GT_Plate4_Multiplex_12_17_2024_Quantification_Cq_Results.csv")
dim (Plate4_MPX)

Plate4_Telo <- read.csv("Telo_mtDNA/Data/GT_Plate4_Telomeres_12_17_2024_Quantification_Cq_Results.csv")
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
write.csv(file = "Telo_mtDNA/Data/Plate4_FinalData.csv", Plate4_FinalData, row.names = FALSE)

############################################
## ## Plate 5##

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate5_MPX <- read.csv("Telo_mtDNA//Data/GT_Plate5_Multiplex_12_17_2024_Quantification_Cq_Results.csv")
dim (Plate5_MPX)

Plate5_Telo <- read.csv("Telo_mtDNA/Data/GT_Plate5_Telomeres_12_17_2024_Quantification_Cq_Results.csv")
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
write.csv(file = "Telo_mtDNA/Data/Plate5_FinalData.csv", Plate5_FinalData, row.names = FALSE)

#######################################
## ## Plate 6##

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate6_MPX <- read.csv("Telo_mtDNA/Data/GT_Plate6_Multiplex_12_18_2024_Quantification_Cq_Results.csv")
dim (Plate6_MPX)

Plate6_MPX <- Plate6_MPX[, -c(1)]

Plate6_Telo <- read.csv("Telo_mtDNA/Data/GT_Plate6_Telomeres_12_18_2024_Quantification_Cq_Results.csv")
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
write.csv(file = "Telo_mtDNA/Data/Plate6_FinalData.csv", Plate6_FinalData, row.names = FALSE)

# Merge all FinalData together, for each plate, by concatenating the rows

# Load each plate's final data
Plate1_FinalData <- read.csv("Telo_mtDNA/Data/Plate1_FinalData.csv")
Plate2_FinalData <- read.csv("Telo_mtDNA/Data/Plate2_FinalData.csv")
Plate3_FinalData <- read.csv("Telo_mtDNA/Data/Plate3_FinalData.csv")
Plate4_FinalData <- read.csv("Telo_mtDNA/Data/Plate4_FinalData.csv")
Plate5_FinalData <- read.csv("Telo_mtDNA/Data/Plate5_FinalData.csv")
Plate6_FinalData <- read.csv("Telo_mtDNA/Data/Plate6_FinalData.csv")

# Concatenate the datasets
qPCR_FinalData <- rbind(Plate1_FinalData, Plate2_FinalData, Plate3_FinalData, Plate4_FinalData, Plate5_FinalData, Plate6_FinalData)

# Verify the dimensions of the combined data
print(dim(qPCR_FinalData))

# Merge Final Data with Trait MetaData for your individuals
# Import .csv files for each run for a particular plate of samples 
Trait <- read.csv("Telo_mtDNA/Data/Trait_MetaData.csv")
dim (Trait)
dim(qPCR_FinalData)

# Merge both datasets
FinalData <- merge(qPCR_FinalData, Trait, by = c("Sample"))

# Save the merged data to a final CSV file
write.csv(FinalData, "Telo_mtDNA/Data/GT_FinalData.csv", row.names = FALSE)

# Optional: Print the dimensions of the final merged data
print(dim(FinalData))

############################################
################ DATA PREPARATION ############
############################################


#################### Calculating the change in Telomere length and mtDNA copy number ####################
# Load the data
data <- read_csv("Telo_mtDNA/Data/GT_FinalData.csv")

# Creating a new column called Tortoise_ID by extracting the individual ID from the column Sample
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

# Calculate changes between timepoints, i.e. absolute difference
## The value at the later timepoint minus the value at the earlier timepoint, using lag() within each tortoise

interval_A_to_B <- interval_A_to_B %>%
  group_by(Tortoise_ID) %>%
  mutate(
    mtDNA_change = mtDNA.Mean - lag(mtDNA.Mean),
    Telomeres_change = Telomeres.per.cell - lag(Telomeres.per.cell),
    Transition = paste0(lag(Time_Point), "_to_", Time_Point)
  ) %>%
  filter(!is.na(mtDNA_change), !is.na(Telomeres_change)) %>%
  ungroup()

# Calculate changes from one timepoint to the next, i.e. from timepoint B to C
interval_B_to_C <- interval_B_to_C %>%
  group_by(Tortoise_ID) %>%
  mutate(
    mtDNA_change = mtDNA.Mean - lag(mtDNA.Mean),
    Telomeres_change = Telomeres.per.cell - lag(Telomeres.per.cell),
    Transition = paste0(lag(Time_Point), "_to_", Time_Point)
  ) %>%
  filter(!is.na(mtDNA_change), !is.na(Telomeres_change)) %>%
  ungroup()

# Calculate changes from one timepoint to the next, i.e. from timepoint C to D
interval_C_to_D <- interval_C_to_D %>%
  group_by(Tortoise_ID) %>%
  mutate(
    mtDNA_change = mtDNA.Mean - lag(mtDNA.Mean),
    Telomeres_change = Telomeres.per.cell - lag(Telomeres.per.cell),
    Transition = paste0(lag(Time_Point), "_to_", Time_Point)
  ) %>%
  filter(!is.na(mtDNA_change), !is.na(Telomeres_change)) %>%
  ungroup()

# Merge all intervals together
merged_intervals <- rbind(interval_A_to_B, interval_B_to_C, interval_C_to_D)

# Now pivot to wide format
change_wide <- merged_intervals %>%
  select(Tortoise_ID, Transition, mtDNA_change, Telomeres_change) %>%
  pivot_wider(
    names_from = Transition,
    values_from = c(mtDNA_change, Telomeres_change),
    names_glue = "{.value}_{Transition}"
  )

# Produce a .csv file with one row per Tortoise_ID and each column represents a change in either mtDNA or telomere length from one timepoint to the next
write_csv(change_wide, "Telo_mtDNA/Data/GT_Change.csv")


###### Joining the Change file with the Growth Rate file to run the mixed-effect model ######

# Removing Not_viable individuals or those with deviations from normal growth
growth_data <- read_csv("Telo_mtDNA/Data/GT_GrowthData.csv") %>%
  filter(!(Tortoise_ID %in% c("Not_Viable", "GT2023_N06.01", "GT2023_N05.06", "GT2023_N15.04"))) %>%
  drop_na()  # removes any remaining rows with NA

# Split the ID into two parts, pad both to two digits, and paste back together
change_wide <- read_csv("Telo_mtDNA/Data/GT_Change.csv") %>%
  mutate(
    Tortoise_ID = sprintf("GT2023_N%02d.%02d",
                          as.integer(sub("\\..*", "", Tortoise_ID)),  # part before the decimal
                          as.integer(sub(".*\\.", "", Tortoise_ID))   # part after the decimal
    )
  )

# Now merge the two datasets by Tortoise_ID
merged_data <- left_join(growth_data, change_wide, by = "Tortoise_ID")

# Produce a merged .csv file containing growth rate and change in mtDNA and telomere length from one timepoint to the next
write_csv(merged_data, "Telo_mtDNA/Data/GT_Growth_Telo_mtDNA.csv")


############################################
### TELOMERES ###
############################################

# Load data
datum <- read.csv("Telo_mtDNA/Data/GT_Growth_Telo_mtDNA.csv", na.strings = "na")

############################################
### 1. PREPARE LONG FORMAT DATA ###
############################################

# Create long format for telomere changes
telomere_long <- datum %>%
  select(Tortoise_ID, Treatment, Nest_ID, Tank,
         Telomeres_change_A_to_B, Telomeres_change_B_to_C, Telomeres_change_C_to_D,
         Growth_rate_During, Growth_rate_3_Weeks_Post, Growth_rate_3_Months_Post) %>%
  pivot_longer(
    cols = c(Telomeres_change_A_to_B, Telomeres_change_B_to_C, Telomeres_change_C_to_D),
    names_to = "Time_Interval",
    values_to = "Telomere_Change"
  ) %>%
  mutate(
    Time_Interval = factor(Time_Interval,
                           levels = c("Telomeres_change_A_to_B", 
                                      "Telomeres_change_B_to_C", 
                                      "Telomeres_change_C_to_D"),
                           labels = c("Winter_2023", "Early_Spring_2024", "Late_Spring_2024")),
    Telomere_Change = as.numeric(Telomere_Change)
  ) %>%
  filter(!is.na(Telomere_Change))

# Add corresponding growth rates
telomere_long <- telomere_long %>%
  mutate(Growth_Rate = case_when(
    Time_Interval == "Winter_2023" ~ Growth_rate_During,
    Time_Interval == "Early_Spring_2024" ~ Growth_rate_3_Weeks_Post,
    Time_Interval == "Late_Spring_2024" ~ Growth_rate_3_Months_Post
  )) %>%
  filter(!is.na(Growth_Rate))

# Check for outliers (|z| > 3) within each interval
telomere_long <- telomere_long %>%
  group_by(Time_Interval) %>%
  mutate(
    z_telo = abs((Telomere_Change - mean(Telomere_Change, na.rm = TRUE)) / 
                   sd(Telomere_Change, na.rm = TRUE)),
    z_growth = abs((Growth_Rate - mean(Growth_Rate, na.rm = TRUE)) / 
                     sd(Growth_Rate, na.rm = TRUE)),
    outlier = z_telo > 3 | z_growth > 3
  ) %>%
  ungroup()

# Print outliers
outliers <- telomere_long %>%
  filter(outlier) %>%
  select(Tortoise_ID, Time_Interval, Telomere_Change, Growth_Rate, z_telo, z_growth)

print(outliers)

## Print sample size by treatment and time interval
telomere_long %>%
  group_by(Time_Interval, Treatment) %>%
  summarise(n = n(), .groups = "drop") %>%
  print()

############################################################################
## ANALYSIS 1: Is the change in telomere length different between treatments?
############################################################################

# Linear model
model_treatment <- lmer(Telomere_Change ~ Treatment * Time_Interval + (1 | Tortoise_ID) + (1 | Nest_ID),
                      data = telomere_long, control = lmerControl(optimizer = "bobyqa"))

summary(model_treatment)
anova(model_treatment)

# Get estimated marginal means for Treatment × Time_Interval
emm <- emmeans(model_treatment, ~ Treatment | Time_Interval)

# Pairwise comparisons within each time interval
pairs(emm, adjust = "tukey")

# Get contrasts for each treatment separately at each interval
emm2 <- emmeans(model_treatment, ~ Time_Interval | Treatment)
pairs(emm2, adjust = "tukey")

confint(pairs(emm, adjust = "tukey"))

### The 95% C.I. for Winter 2023 is ± 4800386 (1995007-(-2805379))
### The 95% C.I. for Early Spring 2024 is ± 4156489 (2970889-(-1185600))
### The 95% C.I. for Late Spring 2024 is ± 4189023 (1100648-(-3088375))

########################################
# Model diagnostics
#######################################

# 1. NORMALITY OF RESIDUALS
par(mfrow = c(1, 2))
qqnorm(resid(model_treatment), main = "Q-Q Plot: Residuals")
qqline(resid(model_treatment), col = "red")
hist(resid(model_treatment), main = "Histogram of Residuals", 
     xlab = "Residuals", breaks = 20)
par(mfrow = c(1, 1))

# Shapiro-Wilk test
shapiro.test(resid(model_treatment))

# 2. HOMOGENEITY OF VARIANCE (CORRECTED)
par(mfrow = c(1, 2))

# Residuals vs Fitted
plot(fitted(model_treatment), resid(model_treatment),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lwd = 2)

# Residuals by Treatment (use boxplot for factors)
boxplot(resid(model_treatment) ~ telomere_long$Treatment,
        xlab = "Treatment", ylab = "Residuals",
        main = "Residuals by Treatment",
        col = c("lightblue", "lightcoral"))
abline(h = 0, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

# Levene's test for homogeneity of variance
leveneTest(resid(model_treatment) ~ Treatment, data = telomere_long)
leveneTest(resid(model_treatment) ~ Time_Interval, data = telomere_long)
leveneTest(resid(model_treatment) ~ Treatment * Time_Interval, data = telomere_long)


###############################################
############## DATA VISUALISATION ##########
##############################################

# Define color-blind friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

######## Plot Change in Telomere Length ##############
# Prepare data for plotting

all_tel_data <- telomere_long %>%
  mutate(
    Time_Point = factor(Time_Interval,
                        levels = c("Winter_2023", "Early_Spring_2024", "Late_Spring_2024"),
                        labels = c("Winter 2023\n(ED – BD)", 
                                   "Early Spring 2024\n(3wkPD – ED)", 
                                   "Late Spring 2024\n(3moPD – 3wkPD)")),
    Treatment = factor(Treatment,
                       levels = c("Cold-Dormancy", "Constant-Warmth"))
  ) %>%
  filter(
    !is.na(Telomere_Change),
    !(sub("_.*", "", Tortoise_ID) %in% c("GT2023_N06.01", "GT2023_N05.06", "GT2023_N15.04"))
  )

# Plot
tel_plot <- ggplot(all_tel_data, aes(x = Time_Point, y = Telomere_Change, color = Treatment)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7,
    position = position_dodge(0.8),
    fill = "white",
    alpha = 0.6,
    linewidth = 0.6
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
    size = 2.5,
    alpha = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.8,
    aes(group = Treatment),
    position = position_dodge(width = 0.8)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20)),
    axis.text.x = element_text(size = 12, hjust = 0.5, lineheight = 1.2),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  scale_color_manual(
    values = c("Cold-Dormancy" = "#0072B2", "Constant-Warmth" = "#D55E00"),
    name = "Treatment",
    labels = c("Cold-Dormancy", "Constant-Warmth")
  ) +
  labs(
    x = "Time Interval",
    y = "Change in telomere length (T/S)"
  )

# Show plot
tel_plot

ggsave(tel_plot, file = "Telo_mtDNA/Figures/Telomeres_Change_Treatment.png", width = 9, height = 7, dpi = 600)

###############################################################################
## ANALYSIS 2: Is the change in telomere length correlated with fast growth rate?
###############################################################################

# Correlation by Treatment
cor_by_treatment <- telomere_long %>%
  group_by(Treatment) %>%
  summarise(
    cor_test = list(cor.test(Telomere_Change, Growth_Rate, method = "spearman", exact = FALSE))
  ) %>%
  mutate(
    rho = sapply(cor_test, \(x) x$estimate),
    p = sapply(cor_test, \(x) x$p.value)
  ) %>%
  select(Treatment, rho, p)

print(cor_by_treatment)

# Correlation by Treatment and Time_Interval
cor_by_treatment_interval <- telomere_long %>%
  group_by(Treatment, Time_Interval) %>%
  summarise(
    n = n(),
    cor_test = list(cor.test(Telomere_Change, Growth_Rate, method = "spearman", exact = FALSE))
  ) %>%
  mutate(
    rho = sapply(cor_test, \(x) x$estimate),
    p = sapply(cor_test, \(x) x$p.value)
  ) %>%
  select(Treatment, Time_Interval, n, rho, p)

print(cor_by_treatment_interval)

###############################################
## VISUALIZATION
###############################################

# Color palette
cbbPalette <- c("Cold-Dormancy" = "#0072B2", "Constant-Warmth" = "#D55E00")

# Determine global Y-axis limits and breaks
all_telomeres <- telomere_long$Telomere_Change

y_step <- 5e6
y_min <- floor(min(all_telomeres, na.rm = TRUE) / y_step) * y_step
y_max <- ceiling(max(all_telomeres, na.rm = TRUE) / y_step) * y_step
y_breaks <- seq(y_min, y_max, by = y_step)

# Standardize formatting
point_size <- 3
point_alpha <- 0.6
axis_title_size <- 18
legend_title_size <- 15
legend_text_size <- 12
line_width <- 1

##############################
# Winter 2023 (A -> B)
##############################
data_winter <- telomere_long %>%
  filter(Time_Interval == "Winter_2023")

Telo_Cor_Winter_2023 <- ggplot(data_winter,
                               aes(x = Growth_Rate,
                                   y = Telomere_Change,
                                   color = Treatment)) +
  geom_point(size = point_size, alpha = point_alpha) +
  geom_smooth(method = "lm", se = FALSE, linewidth = line_width) +
  scale_color_manual(values = cbbPalette) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
  scale_y_continuous(limits = c(y_min, y_max),
                     breaks = y_breaks,
                     labels = label_scientific()) +
  labs(x = "Growth rate (g/day)\nWinter 2023 (ED-BD)",
       y = "Change in telomere length (T/S)") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = axis_title_size, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = axis_title_size, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = legend_title_size, face = "bold"),
    legend.text = element_text(size = legend_text_size),
    legend.position = c(0.85, 0.95),  # Inside top-right
    legend.justification = c(0.5, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  )
##############################
# Early Spring 2024 (B -> C)
##############################
data_early_spring <- telomere_long %>%
  filter(Time_Interval == "Early_Spring_2024")

Telo_Cor_Early_Spring_2024 <- ggplot(data_early_spring,
                                     aes(x = Growth_Rate,
                                         y = Telomere_Change,
                                         color = Treatment)) +
  geom_point(size = point_size, alpha = point_alpha) +
  geom_smooth(method = "lm", se = FALSE, linewidth = line_width) +
  scale_color_manual(values = cbbPalette) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
  scale_y_continuous(limits = c(y_min, y_max),
                     breaks = y_breaks,
                     labels = label_scientific()) +
  labs(x = "Growth rate (g/day)\nEarly Spring 2024 (3wkPD – ED)",
       y = "Change in telomere length (T/S)") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = axis_title_size, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = axis_title_size, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = legend_title_size, face = "bold"),
    legend.text = element_text(size = legend_text_size),
    legend.position = "none",
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  )

##############################
# Late Spring 2024 (C -> D)
##############################
data_late_spring <- telomere_long %>%
  filter(Time_Interval == "Late_Spring_2024")

Telo_Cor_Late_Spring_2024 <- ggplot(data_late_spring,
                                    aes(x = Growth_Rate,
                                        y = Telomere_Change,
                                        color = Treatment)) +
  geom_point(size = point_size, alpha = point_alpha) +
  geom_smooth(method = "lm", se = FALSE, linewidth = line_width) +
  scale_color_manual(values = cbbPalette) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
  scale_y_continuous(limits = c(y_min, y_max),
                     breaks = y_breaks,
                     labels = label_scientific()) +
  labs(x = "Growth rate (g/day)\nLate Spring 2024 (3moPD – 3wkPD)",
       y = "Change in telomere length (T/S)") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = axis_title_size, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = axis_title_size, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = legend_title_size, face = "bold"),
    legend.text = element_text(size = legend_text_size),
    legend.position = "none",
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  )

library(patchwork)

# Print plots
Telo_Cor_Winter_2023
Telo_Cor_Early_Spring_2024
Telo_Cor_Late_Spring_2024

# Combine plots with bold tags inside panels (top-left)
combined_plot <- (Telo_Cor_Winter_2023 + theme(plot.tag.position = c(0.2, 1))) / 
  (Telo_Cor_Early_Spring_2024 + theme(plot.tag.position = c(0.2, 1))) / 
  (Telo_Cor_Late_Spring_2024 + theme(plot.tag.position = c(0.2, 1))) +
  plot_annotation(
    tag_levels = 'A'
  ) &
  theme(plot.tag = element_text(face = "bold", size = 20))

combined_plot

# Save combined plot
ggsave(combined_plot, file = "Telo_mtDNA/Figures/Telomeres_Correlation_Growth_Rate.png", width = 9, height = 15, dpi = 600)

#########################################
### Mitochondrial DNA ###
#########################################

# Read and prepare data
mtdna_data <- read_csv("Telo_mtDNA//Data/GT_Merged_Clean.csv") %>%
  filter(!(Tortoise_ID %in% c("GT2023_N06.01", "GT2023_N05.06", "GT2023_N15.04"))) %>%
  mutate(
    Time_Point = factor(Time_Point, 
                        levels = c("A", "B", "C", "D"),
                        labels = c("Before", "Dormancy_end", "3-Wk.Post", "3-Mo.Post")),
    Treatment = factor(Treatment, 
                       levels = c("Cold-Dormancy", "Constant-Warmth")),
    Nest_ID = factor(Nest_ID),
    Tank = factor(Tank),
    Tortoise_ID = factor(Tortoise_ID),
    Plate_ID = factor(Plate_ID)
  )

# Fit model
model_mtdna <- lmer(mtDNA.Mean ~ Treatment * Time_Point + (1 | Tortoise_ID) + (1 | Plate_ID),
                    data = mtdna_data,
                    control = lmerControl(optimizer = "bobyqa"))

str(mtdna_data)

summary(model_mtdna)
anova(model_mtdna)

## Post-hoc comparisons ##
emmeans_mtdna <- emmeans(model_mtdna, pairwise ~ Treatment | Time_Point)

# Compare timepoints within each treatment
emmeans_time <- emmeans(model_mtdna, pairwise ~ Time_Point | Treatment)
summary(emmeans_time)

# Compare between treatments at each timepoint
emmeans_treatment_by_time <- emmeans(model_mtdna, ~ Treatment | Time_Point)
pairs(emmeans_treatment_by_time)

confint(pairs(emmeans_treatment_by_time))

### The 95% C.I. for End-Dormancy is ± 4.57 (-0.45-(-5.020))
### The 95% C.I. for 3-Week Post-Dormancy is ± 4.57 (5.49-(0.916))
### The 95% C.I. for 3-Month Post-Dormancy ± 4.67 (2.70-(-1.970))

################################################
############## Model diagnostics ###############
################################################

# 1. NORMALITY OF RESIDUALS
par(mfrow = c(1, 2))
qqnorm(resid(model_mtdna), main = "Q-Q Plot: Residuals")
qqline(resid(model_mtdna), col = "red")
hist(resid(model_mtdna), main = "Histogram of Residuals", 
     xlab = "Residuals", breaks = 20)
par(mfrow = c(1, 1))

shapiro.test(resid(model_mtdna))

# 2. HOMOGENEITY OF VARIANCE
par(mfrow = c(1, 2))

plot(fitted(model_mtdna), resid(model_mtdna),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lwd = 2)

boxplot(resid(model_mtdna) ~ mtdna_data$Treatment,
        xlab = "Treatment", ylab = "Residuals",
        main = "Residuals by Treatment",
        col = c("#0072B2", "#D55E00"))
abline(h = 0, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

leveneTest(resid(model_mtdna) ~ Treatment, data = mtdna_data)
leveneTest(resid(model_mtdna) ~ Time_Point, data = mtdna_data)

# 3. CHECK FOR OUTLIERS
mtdna_data <- mtdna_data %>%
  group_by(Time_Point) %>%
  mutate(
    z_mtdna = abs((mtDNA.Mean - mean(mtDNA.Mean, na.rm = TRUE)) / 
                    sd(mtDNA.Mean, na.rm = TRUE)),
    outlier = z_mtdna > 3
  ) %>%
  ungroup()

# Print outliers
outliers_mtdna <- mtdna_data %>%
  filter(outlier) %>%
  select(Sample, Tortoise_ID, Time_Point, Treatment, mtDNA.Mean, z_mtdna)

print(outliers_mtdna)

# Summary of outliers by time point
mtdna_data %>%
  group_by(Time_Point) %>%
  summarise(
    n_total = n(),
    n_outliers = sum(outlier, na.rm = TRUE),
    percent_outliers = round(100 * n_outliers / n_total, 1)
  )


###############################################
############## DATA VISUALISATION ##########
##############################################

# Define color-blind friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

######## Plot mtDNA Copy Number ##############

# Prepare data for plotting
all_mtdna_data <- mtdna_data %>%
  mutate(
    Time_Point = factor(Time_Point,
                        levels = c("Before", "Dormancy_end", "3-Wk.Post", "3-Mo.Post"),
                        labels = c("BD", "ED", "3wkPD", "3moPD")),
    Treatment = factor(Treatment,
                       levels = c("Cold-Dormancy", "Constant-Warmth"))
  ) %>%
  filter(!is.na(mtDNA.Mean))

## Plot ##

# Define y positions for significance brackets
ypos <- max(all_mtdna_data$mtDNA.Mean, na.rm = TRUE) + 3

mtdna_plot <- ggplot() +
  # BD - single gray box (combined treatments)
  geom_boxplot(data = subset(all_mtdna_data, Time_Point == "BD"),
               aes(x = Time_Point, y = mtDNA.Mean),
               color = "black", fill = "gray70", outlier.shape = NA, 
               width = 0.6, linewidth = 0.6) +
  geom_point(data = subset(all_mtdna_data, Time_Point == "BD"),
             aes(x = Time_Point, y = mtDNA.Mean),
             position = position_jitter(width = 0.1), 
             size = 2, color = "black", fill = "white", shape = 21, stroke = 0.5) +
  
  # Mean point for BD (black, single combined point)
  stat_summary(data = subset(all_mtdna_data, Time_Point == "BD"),
               aes(x = Time_Point, y = mtDNA.Mean),
               fun = mean, geom = "point", size = 2, color = "black") +
  
  # Other timepoints - separate boxes by Treatment (white fill, colored outlines)
  geom_boxplot(data = subset(all_mtdna_data, Time_Point != "BD"),
               aes(x = Time_Point, y = mtDNA.Mean, color = Treatment),
               fill = "white", outlier.shape = NA, width = 0.6, 
               position = position_dodge(0.8), linewidth = 0.6) +
  geom_jitter(data = subset(all_mtdna_data, Time_Point != "BD"),
              aes(x = Time_Point, y = mtDNA.Mean, color = Treatment),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
              size = 2, alpha = 0.4) +
  
  # Mean points for other timepoints
  stat_summary(data = subset(all_mtdna_data, Time_Point != "BD"),
               aes(x = Time_Point, y = mtDNA.Mean, color = Treatment, group = Treatment),
               fun = mean, geom = "point", size = 2,
               position = position_dodge(width = 0.8)) +
  
  # Between-treatment significance bars (BLACK horizontal bars)
  # ED: * (p = 0.0192)
  annotate("segment", x = 1.9, xend = 2.1, y = ypos - 3, yend = ypos - 3,
           color = "black", linewidth = 0.8) +
  annotate("text", x = 2.0, y = ypos - 1, label = "*", 
           size = 6, color = "black") +
  
  # 3wkPD: ** (p = 0.0063)
  annotate("segment", x = 2.9, xend = 3.1, y = ypos - 3, yend = ypos - 3,
           color = "black", linewidth = 0.8) +
  annotate("text", x = 3.0, y = ypos - 1, label = "**", 
           size = 6, color = "black") +
  
  # ORANGE BRACKET 1: Constant-Warmth ED → 3wkPD (p < 0.0001)
  annotate("segment", x = 2.2, xend = 2.2, y = ypos, yend = ypos + 3,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("segment", x = 2.2, xend = 3.2, y = ypos + 3, yend = ypos + 3,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("segment", x = 3.2, xend = 3.2, y = ypos, yend = ypos + 3,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("text", x = 2.7, y = ypos + 5, label = "***", 
           size = 6, color = cbbPalette[7]) +
  
  # ORANGE BRACKET 2: Constant-Warmth ED → 3moPD (p = 0.0012)
  annotate("segment", x = 2.2, xend = 2.2, y = ypos + 7, yend = ypos + 10,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("segment", x = 2.2, xend = 4.2, y = ypos + 10, yend = ypos + 10,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("segment", x = 4.2, xend = 4.2, y = ypos + 7, yend = ypos + 10,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("text", x = 3.2, y = ypos + 12, label = "**", 
           size = 6, color = cbbPalette[7]) +
  
  # Theme and labels
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20, b = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "right",
    legend.title.align = 0.5
  ) +
  
  scale_x_discrete(expand = expansion(add = 0.8)) +
  
  # Colors matching your colorblind palette
  scale_color_manual(
    values = c("Cold-Dormancy" = cbbPalette[6], "Constant-Warmth" = cbbPalette[7]),
    name = "Treatment",
    labels = c("Cold-Dormancy", "Constant-Warmth")
  ) +
  
  labs(
    x = "Timepoint",
    y = "mtDNA density (mtDNA/S)"
  )

mtdna_plot

# Save the plot
ggsave(mtdna_plot, file = "Telo_mtDNA/Figures/mtDNA_Treatment.png", width = 9, height = 7, dpi = 600)

# Sample size by Treatment and Time_Point
sample_sizes <- mtdna_data %>%
  group_by(Treatment, Time_Point) %>%
  summarise(n = n(), .groups = 'drop')

print(sample_sizes)

###### Saving the mtDNA_Treatment.png as an object so I can reload it in the Metabolites script and
# combine the metabolites figures into one 4-panel plot

saveRDS(mtdna_plot, file = "Telo_mtDNA/Figures/mtdna_plot.rds")


###########################################################
######### mtDNA change in response to treatment across intervals #########
###########################################################

### 1. PREPARE LONG FORMAT DATA FOR mtDNA CHANGE ###
############################################

# Create long format for mtDNA changes
mtdna_change_long <- datum %>%
  select(Tortoise_ID, Treatment, Nest_ID, Tank,
         mtDNA_change_A_to_B, mtDNA_change_B_to_C, mtDNA_change_C_to_D,
         Growth_rate_During, Growth_rate_3_Weeks_Post, Growth_rate_3_Months_Post) %>%
  pivot_longer(
    cols = c(mtDNA_change_A_to_B, mtDNA_change_B_to_C, mtDNA_change_C_to_D),
    names_to = "Time_Interval",
    values_to = "mtDNA_Change"
  ) %>%
  mutate(
    Time_Interval = factor(Time_Interval,
                           levels = c("mtDNA_change_A_to_B", 
                                      "mtDNA_change_B_to_C", 
                                      "mtDNA_change_C_to_D"),
                           labels = c("Winter_2023", "Early_Spring_2024", "Late_Spring_2024")),
    mtDNA_Change = as.numeric(mtDNA_Change)
  ) %>%
  filter(!is.na(mtDNA_Change))

# Add corresponding growth rates
mtdna_change_long <- mtdna_change_long %>%
  mutate(Growth_Rate = case_when(
    Time_Interval == "Winter_2023" ~ Growth_rate_During,
    Time_Interval == "Early_Spring_2024" ~ Growth_rate_3_Weeks_Post,
    Time_Interval == "Late_Spring_2024" ~ Growth_rate_3_Months_Post
  )) %>%
  filter(!is.na(Growth_Rate))

# Check for outliers (|z| > 3) within each interval
mtdna_change_long <- mtdna_change_long %>%
  group_by(Time_Interval) %>%
  mutate(
    z_mtdna = abs((mtDNA_Change - mean(mtDNA_Change, na.rm = TRUE)) / 
                    sd(mtDNA_Change, na.rm = TRUE)),
    z_growth = abs((Growth_Rate - mean(Growth_Rate, na.rm = TRUE)) / 
                     sd(Growth_Rate, na.rm = TRUE)),
    outlier = z_mtdna > 3 | z_growth > 3
  ) %>%
  ungroup()

# Print outliers
outliers_mtdna <- mtdna_change_long %>%
  filter(outlier) %>%
  select(Tortoise_ID, Time_Interval, mtDNA_Change, Growth_Rate, z_mtdna, z_growth)
print(outliers_mtdna)

# Print sample size by treatment and time interval
mtdna_change_long %>%
  group_by(Time_Interval, Treatment) %>%
  summarise(n = n(), .groups = "drop") %>%
  print()

############################################################################
## ANALYSIS 2: Is the change in mtDNA copy number different between treatments?
############################################################################

# Linear model
model_mtdna_change <- lmer(mtDNA_Change ~ Treatment * Time_Interval + (1 | Tortoise_ID) + (1 | Nest_ID),
                           data = mtdna_change_long, control = lmerControl(optimizer = "bobyqa"))
summary(model_mtdna_change)
anova(model_mtdna_change)

# Get estimated marginal means for Treatment × Time_Interval
emm_mtdna <- emmeans(model_mtdna_change, ~ Treatment | Time_Interval)

# Pairwise comparisons within each time interval
pairs(emm_mtdna, adjust = "tukey")

# Get contrasts for each treatment separately at each interval
emm2_mtdna <- emmeans(model_mtdna_change, ~ Time_Interval | Treatment)
pairs(emm2_mtdna, adjust = "tukey")

########################################
# Model diagnostics
#######################################

# 1. NORMALITY OF RESIDUALS
par(mfrow = c(1, 2))
qqnorm(resid(model_mtdna_change), main = "Q-Q Plot: Residuals")
qqline(resid(model_mtdna_change), col = "red")
hist(resid(model_mtdna_change), main = "Histogram of Residuals", 
     xlab = "Residuals", breaks = 20)
par(mfrow = c(1, 1))

# Shapiro-Wilk test
shapiro.test(resid(model_mtdna_change))

# 2. HOMOGENEITY OF VARIANCE
par(mfrow = c(1, 2))
# Residuals vs Fitted
plot(fitted(model_mtdna_change), resid(model_mtdna_change),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lwd = 2)

# Residuals by Treatment
boxplot(resid(model_mtdna_change) ~ mtdna_change_long$Treatment,
        xlab = "Treatment", ylab = "Residuals",
        main = "Residuals by Treatment",
        col = c("lightblue", "lightcoral"))
abline(h = 0, col = "red", lwd = 2, lty = 2)
par(mfrow = c(1, 1))

# Levene's test for homogeneity of variance
leveneTest(resid(model_mtdna_change) ~ Treatment, data = mtdna_change_long)
leveneTest(resid(model_mtdna_change) ~ Time_Interval, data = mtdna_change_long)
leveneTest(resid(model_mtdna_change) ~ Treatment * Time_Interval, data = mtdna_change_long)


###############################################
############## DATA VISUALISATION ##########
##############################################


# Prepare data for plotting
all_mtdna_change_data <- mtdna_change_long %>%
  mutate(
    Time_Point = factor(Time_Interval,
                        levels = c("Winter_2023", "Early_Spring_2024", "Late_Spring_2024"),
                        labels = c("Winter 2023\n(ED – BD)", 
                                   "Early Spring 2024\n(3wkPD – ED)", 
                                   "Late Spring 2024\n(3moPD – 3wkPD)")),
    Treatment = factor(Treatment,
                       levels = c("Cold-Dormancy", "Constant-Warmth"))
  ) %>%
  filter(
    !is.na(mtDNA_Change),
    !(sub("_.*", "", Tortoise_ID) %in% c("GT2023_N06.01", "GT2023_N05.06", "GT2023_N15.04"))
  )

# Define y positions for significance brackets
ypos <- max(all_mtdna_change_data$mtDNA_Change, na.rm = TRUE) + 3

# Plot
mtdna_change_plot <- ggplot(all_mtdna_change_data, aes(x = Time_Point, y = mtDNA_Change, color = Treatment)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7,
    position = position_dodge(0.8),
    fill = "white",
    alpha = 0.6,
    linewidth = 0.6
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
    size = 2.5,
    alpha = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.8,
    aes(group = Treatment),
    position = position_dodge(width = 0.8)
  ) +
  
  # Between-treatment significance bars (BLACK horizontal bars)
  # Winter 2023: * (p = 0.0245)
  annotate("segment", x = 0.9, xend = 1.1, y = ypos - 3, yend = ypos - 3,
           color = "black", linewidth = 0.8) +
  annotate("text", x = 1.0, y = ypos - 1, label = "*",
           size = 6, color = "black") +
  
  # Early Spring 2024: *** (p = 0.0008)
  annotate("segment", x = 1.9, xend = 2.1, y = ypos - 3, yend = ypos - 3,
           color = "black", linewidth = 0.8) +
  annotate("text", x = 2.0, y = ypos - 1, label = "***",
           size = 6, color = "black") +
  
  # ORANGE BRACKET 1: Constant-Warmth Winter 2023 → Early Spring 2024 (p < 0.0001)
  annotate("segment", x = 1.2, xend = 1.2, y = ypos, yend = ypos + 3,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("segment", x = 1.2, xend = 2.2, y = ypos + 3, yend = ypos + 3,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("segment", x = 2.2, xend = 2.2, y = ypos, yend = ypos + 3,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("text", x = 1.7, y = ypos + 5, label = "***",
           size = 6, color = cbbPalette[7]) +
  
  # ORANGE BRACKET 2: Constant-Warmth Early Spring 2024 → Late Spring 2024 (p = 0.0011)
  annotate("segment", x = 2.2, xend = 2.2, y = ypos + 7, yend = ypos + 10,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("segment", x = 2.2, xend = 3.2, y = ypos + 10, yend = ypos + 10,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("segment", x = 3.2, xend = 3.2, y = ypos + 7, yend = ypos + 10,
           color = cbbPalette[7], linewidth = 0.8) +
  annotate("text", x = 2.7, y = ypos + 12, label = "**",
           size = 6, color = cbbPalette[7]) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20)),
    axis.text.x = element_text(size = 12, hjust = 0.5, lineheight = 1.2),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  scale_color_manual(
    values = c("Cold-Dormancy" = cbbPalette[6], "Constant-Warmth" = cbbPalette[7]),
    name = "Treatment",
    labels = c("Cold-Dormancy", "Constant-Warmth")
  ) +
  labs(
    x = "Time Interval",
    y = "Change in mtDNA density (mtDNA/S)"
  )

mtdna_change_plot

ggsave(mtdna_change_plot, file = "Telo_mtDNA/Figures/mtDNA_Change_Treatment.png",
       width = 9, height = 7, dpi = 600)

###############################################################################
## ANALYSIS 3: Is the change in mtDNA copy number correlated with growth rate?

###############################################################################
# Correlation by Treatment
cor_by_treatment_mtdna <- mtdna_change_long %>%
  group_by(Treatment) %>%
  summarise(
    cor_test = list(cor.test(mtDNA_Change, Growth_Rate, method = "spearman", exact = FALSE))
  ) %>%
  mutate(
    rho = sapply(cor_test, \(x) x$estimate),
    p = sapply(cor_test, \(x) x$p.value)
  ) %>%
  select(Treatment, rho, p)
print(cor_by_treatment_mtdna)

# Correlation by Treatment and Time_Interval
cor_by_treatment_interval_mtdna <- mtdna_change_long %>%
  group_by(Treatment, Time_Interval) %>%
  summarise(
    n = n(),
    cor_test = list(cor.test(mtDNA_Change, Growth_Rate, method = "spearman", exact = FALSE))
  ) %>%
  mutate(
    rho = sapply(cor_test, \(x) x$estimate),
    p = sapply(cor_test, \(x) x$p.value)
  ) %>%
  select(Treatment, Time_Interval, n, rho, p)
print(cor_by_treatment_interval_mtdna)


###############################################
## VISUALIZATION - mtDNA Change vs Growth Rate
###############################################

# Color palette
cbbPalette <- c("Cold-Dormancy" = "#0072B2", "Constant-Warmth" = "#D55E00")

# Determine global Y-axis limits and breaks
all_mtdna <- mtdna_change_long$mtDNA_Change

y_step <- 5
y_min <- floor(min(all_mtdna, na.rm = TRUE) / y_step) * y_step
y_max <- ceiling(max(all_mtdna, na.rm = TRUE) / y_step) * y_step
y_breaks <- seq(y_min, y_max, by = y_step)

# Standardize formatting
point_size <- 3
point_alpha <- 0.6
axis_title_size <- 18
legend_title_size <- 15
legend_text_size <- 12
line_width <- 1

##############################
# Winter 2023 (A -> B)
##############################
data_winter_mtdna <- mtdna_change_long %>%
  filter(Time_Interval == "Winter_2023")

mtDNA_Cor_Winter_2023 <- ggplot(data_winter_mtdna,
                                aes(x = Growth_Rate,
                                    y = mtDNA_Change,
                                    color = Treatment)) +
  geom_point(size = point_size, alpha = point_alpha) +
  geom_smooth(method = "lm", se = FALSE, linewidth = line_width) +
  scale_color_manual(values = cbbPalette) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = y_breaks) +
  labs(x = "Growth rate (g/day)\nWinter 2023 (ED – BD)",
       y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = axis_title_size, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = axis_title_size, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = legend_title_size, face = "bold"),
    legend.text = element_text(size = legend_text_size),
    legend.position = c(0.85, 0.95),
    legend.justification = c(0.5, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  )

##############################
# Early Spring 2024 (B -> C)
##############################
data_early_spring_mtdna <- mtdna_change_long %>%
  filter(Time_Interval == "Early_Spring_2024")

mtDNA_Cor_Early_Spring_2024 <- ggplot(data_early_spring_mtdna,
                                      aes(x = Growth_Rate,
                                          y = mtDNA_Change,
                                          color = Treatment)) +
  geom_point(size = point_size, alpha = point_alpha) +
  geom_smooth(method = "lm", se = FALSE, linewidth = line_width) +
  scale_color_manual(values = cbbPalette) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = y_breaks) +
  labs(x = "Growth rate (g/day)\nEarly Spring 2024 (3wkPD – ED)",
       y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = axis_title_size, face = "bold", margin = margin(t = 20)),
    axis.text = element_text(size = 12),
    legend.position = "none",
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  )

##############################
# Late Spring 2024 (C -> D)
##############################
data_late_spring_mtdna <- mtdna_change_long %>%
  filter(Time_Interval == "Late_Spring_2024")

mtDNA_Cor_Late_Spring_2024 <- ggplot(data_late_spring_mtdna,
                                     aes(x = Growth_Rate,
                                         y = mtDNA_Change,
                                         color = Treatment)) +
  geom_point(size = point_size, alpha = point_alpha) +
  geom_smooth(method = "lm", se = FALSE, linewidth = line_width) +
  scale_color_manual(values = cbbPalette) +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = y_breaks) +
  labs(x = "Growth rate (g/day)\nLate Spring 2024 (3moPD – 3wkPD)",
       y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = axis_title_size, face = "bold", margin = margin(t = 20)),
    axis.text = element_text(size = 12),
    legend.position = "none",
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  )

library(patchwork)

# Combine plots
combined_plot_mtdna <- (mtDNA_Cor_Winter_2023 + theme(plot.tag.position = c(0.2, 1))) /
  (mtDNA_Cor_Early_Spring_2024 + theme(plot.tag.position = c(0.2, 1))) /
  (mtDNA_Cor_Late_Spring_2024 + theme(plot.tag.position = c(0.2, 1))) +
  plot_annotation(
    tag_levels = 'A'
  ) &
  theme(plot.tag = element_text(face = "bold", size = 20))

# Wrap and add shared y-axis label
combined_plot_mtdna <- wrap_elements(combined_plot_mtdna) +
  labs(tag = "Change in mtDNA density (mtDNA/S)") +
  theme(
    plot.tag = element_text(size = axis_title_size, face = "bold", angle = 90),
    plot.tag.position = "left"
  )

combined_plot_mtdna

# Save combined plot
ggsave(combined_plot_mtdna, file = "Telo_mtDNA/Figures/mtDNA_Correlation_Growth_Rate.png",
       width = 9, height = 15, dpi = 600)

######## Combining both Telomeres_Correlation_Growth_Rate and mtDNA_Correlation_Growth_Rate in one panel ###
library(cowplot)
library(grid)

# Helper theme for tag styling
tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 19),
  plot.tag.position = c(-0.02, 0.98)
)

# --- Row 1: Winter 2023 ---
p1_telo <- Telo_Cor_Winter_2023 +
  labs(x = NULL, y = NULL, tag = "A") +
  theme(legend.position = "none",
        plot.margin = margin(5, 30, 5, 10)) +
  tag_theme

p1_mtdna <- mtDNA_Cor_Winter_2023 +
  labs(x = NULL, y = NULL, tag = "D") +
  theme(legend.position = c(0.98, 0.06),
        legend.justification = c(1, 0),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1.2, "lines"),
        plot.margin = margin(5, 10, 5, 30)) +
  tag_theme

# --- Row 2: Early Spring 2024 ---
p2_telo <- Telo_Cor_Early_Spring_2024 +
  labs(x = NULL, y = NULL, tag = "B") +
  theme(legend.position = "none",
        plot.margin = margin(5, 30, 5, 10)) +
  tag_theme

p2_mtdna <- mtDNA_Cor_Early_Spring_2024 +
  labs(x = NULL, y = NULL, tag = "E") +
  theme(legend.position = "none",
        plot.margin = margin(5, 10, 5, 30)) +
  tag_theme

# --- Row 3: Late Spring 2024 ---
p3_telo <- Telo_Cor_Late_Spring_2024 +
  labs(x = NULL, y = NULL, tag = "C") +
  theme(legend.position = "none",
        plot.margin = margin(5, 30, 5, 10)) +
  tag_theme

p3_mtdna <- mtDNA_Cor_Late_Spring_2024 +
  labs(x = NULL, y = NULL, tag = "F") +
  theme(legend.position = "none",
        plot.margin = margin(5, 10, 5, 30)) +
  tag_theme

# Column headers
col_header_telo <- wrap_elements(
  grid::textGrob("Change in telomere length (T/S)",
                 gp = grid::gpar(fontsize = 18, fontface = "bold"))
)

col_header_mtdna <- wrap_elements(
  grid::textGrob("Change in mtDNA density (mtDNA/S)",
                 x = 0.55,
                 gp = grid::gpar(fontsize = 18, fontface = "bold"))
)
# Row labels (rotated, pushed close to panels)
row_label_winter <- wrap_elements(
  grid::textGrob("Winter 2023\n(ED – BD)",
                 rot = 90,
                 x = 0.9,
                 gp = grid::gpar(fontsize = 21, fontface = "bold"))
)

row_label_early <- wrap_elements(
  grid::textGrob("Early Spring 2024\n(3wkPD – ED)",
                 rot = 90,
                 x = 0.9,
                 gp = grid::gpar(fontsize = 21, fontface = "bold"))
)

row_label_late <- wrap_elements(
  grid::textGrob("Late Spring 2024\n(3moPD – 3wkPD)",
                 rot = 90,
                 x = 0.9,
                 gp = grid::gpar(fontsize = 21, fontface = "bold"))
)

# Empty spacer for top-left corner
spacer <- plot_spacer()

# Build the main patchwork grid (without x-axis label)
main_grid <-
  (spacer | col_header_telo | col_header_mtdna) /
  (row_label_winter | p1_telo | p1_mtdna) /
  (row_label_early  | p2_telo | p2_mtdna) /
  (row_label_late   | p3_telo | p3_mtdna) +
  plot_layout(
    heights = c(0.08, 1, 1, 1),
    widths  = c(0.3, 1, 1)  # <-- only change: reduced from 1 to 0.3
  ) &
  theme(plot.background = element_blank(),
        panel.background = element_blank())

# Use cowplot to add a centered x-axis label below the entire grid
x_label <- ggdraw() +
  draw_label("Growth rate (g/day)", 
             x = 0.67, y = 0.5,
             fontface = "bold", size = 21) +
  theme_void() +
  theme(plot.background = element_blank(),
        panel.background = element_blank())

# Stack the main grid with the x-axis label below it
final_plot <- plot_grid(
  main_grid,
  x_label,
  ncol = 1,
  rel_heights = c(1, 0.03)
) + theme(plot.margin = margin(5, 10, 5, -10))

final_plot

# Save
ggsave(final_plot,
       file = "Telo_mtDNA/Figures/Combined_Telomere_mtDNA_Correlation_Growth_Rate.png",
       width = 16, height = 15, dpi = 600,
       bg = "white")