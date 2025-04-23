setwd("/Users/anetfilipova/Library/CloudStorage/Box-Box/TS_Lab_AnoleAging/NIH_R15_Anole_ModelAging_IGFs/Cox-Lab_Adult_KnownAges/AgeXSex_Telomeres/Data/Data_Analysis/Raw_Data_CqValues/Plate1")
library(dplyr)
library(stringr)

# Clear memory
rm(list=ls(all = TRUE))


#Process EDIT all the runs for EACH Plate of samples sepearately
#Use the ...Quantification Cq Results.csv files
#Copy and paste the block of code and repeat on each plate of samples.
#You can use find replace for the dataset name, Plate1 for Plate2


## Plate 1##
# Import .csv files for each run for a particular plate of samples 
Plate1_MPX <- read.csv("CKA_Multiplex_First_24_qPCR_CqResults.csv")
dim (Plate1_MPX)

Plate1_Telo <- read.csv("CKA_Telomeres_First_24_qPCR_CqResults.csv")
dim (Plate1_Telo)

# Concatenate data across runs for the same plate of samples
Plate1<- rbind(Plate1_MPX, Plate1_Telo)
dim(Plate1)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate1$PlateID <- "Plate1"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate1$Target[Plate1$Fluor == "Cy5"] <- "scnag"
Plate1$Target[Plate1$Fluor == "HEX"] <- "mtdna"
Plate1$Target[Plate1$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate1$Diff_AVG_Cq <- abs(Plate1$Cq - Plate1$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.3 threshold
Plate1$Flag_outlier <- ifelse(Plate1$Target %in% c("mtdna", "scnag") & Plate1$Diff_AVG_Cq > 0.3, "yes", "no")

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
  rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)

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

## ## Plate 2##
setwd("/Users/anetfilipova/Library/CloudStorage/Box-Box/TS_Lab_AnoleAging/NIH_R15_Anole_ModelAging_IGFs/Cox-Lab_Adult_KnownAges/AgeXSex_Telomeres/Data/Data_Analysis/Raw_Data_CqValues/Plate2")
library(dplyr)
library(stringr)

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate2_MPX <- read.csv("CKA_Multiplex_Second_24_qPCR_CqResults.csv")
dim (Plate2_MPX)

Plate2_Telo <- read.csv("CKA_Telomeres_Second_24_qPCR_CqResults.csv")
dim (Plate2_Telo)

# Concatenate data across runs for the same plate of samples
Plate2<- rbind(Plate2_MPX, Plate2_Telo)
dim(Plate2)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate2$PlateID <- "Plate2"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate2$Target[Plate2$Fluor == "Cy5"] <- "scnag"
Plate2$Target[Plate2$Fluor == "HEX"] <- "mtdna"
Plate2$Target[Plate2$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate2$Diff_AVG_Cq <- abs(Plate2$Cq - Plate2$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.3 threshold
Plate2$Flag_outlier <- ifelse(Plate2$Target %in% c("mtdna", "scnag") & Plate2$Diff_AVG_Cq > 0.3, "yes", "no")

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
  rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)

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
setwd("/Users/anetfilipova/Library/CloudStorage/Box-Box/TS_Lab_AnoleAging/NIH_R15_Anole_ModelAging_IGFs/Cox-Lab_Adult_KnownAges/AgeXSex_Telomeres/Data/Data_Analysis/Raw_Data_CqValues/Plate3")
library(dplyr)
library(stringr)

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate3_MPX <- read.csv("CKA_Multiplex_Third_24_qPCR_CqResults.csv")
dim (Plate3_MPX)

Plate3_Telo <- read.csv("CKA_Telomeres_Third_24_qPCR_CqResults.csv")
dim (Plate3_Telo)

# Concatenate data across runs for the same plate of samples
Plate3<- rbind(Plate3_MPX, Plate3_Telo)
dim(Plate3)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate3$PlateID <- "Plate3"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate3$Target[Plate3$Fluor == "Cy5"] <- "scnag"
Plate3$Target[Plate3$Fluor == "HEX"] <- "mtdna"
Plate3$Target[Plate3$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate3$Diff_AVG_Cq <- abs(Plate3$Cq - Plate3$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.3 threshold
Plate3$Flag_outlier <- ifelse(Plate3$Target %in% c("mtdna", "scnag") & Plate3$Diff_AVG_Cq > 0.3, "yes", "no")

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
  rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)

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
setwd("/Users/anetfilipova/Library/CloudStorage/Box-Box/TS_Lab_AnoleAging/NIH_R15_Anole_ModelAging_IGFs/Cox-Lab_Adult_KnownAges/AgeXSex_Telomeres/Data/Data_Analysis/Raw_Data_CqValues/Plate4")
library(dplyr)
library(stringr)

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate4_MPX <- read.csv("CKA_Multiplex_Fourth_24_qPCR_CqResults.csv")
dim (Plate4_MPX)

Plate4_Telo <- read.csv("CKA_Telomeres_Fourth_24_qPCR_CqResults.csv")
dim (Plate4_Telo)

# Concatenate data across runs for the same plate of samples
Plate4<- rbind(Plate4_MPX, Plate4_Telo)
dim(Plate4)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate4$PlateID <- "Plate4"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate4$Target[Plate4$Fluor == "Cy5"] <- "scnag"
Plate4$Target[Plate4$Fluor == "HEX"] <- "mtdna"
Plate4$Target[Plate4$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate4$Diff_AVG_Cq <- abs(Plate4$Cq - Plate4$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.3 threshold
Plate4$Flag_outlier <- ifelse(Plate4$Target %in% c("mtdna", "scnag") & Plate4$Diff_AVG_Cq > 0.3, "yes", "no")

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
  rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)

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

## ## Plate 5##
setwd("/Users/anetfilipova/Library/CloudStorage/Box-Box/TS_Lab_AnoleAging/NIH_R15_Anole_ModelAging_IGFs/Cox-Lab_Adult_KnownAges/AgeXSex_Telomeres/Data/Data_Analysis/Raw_Data_CqValues/Plate5")
library(dplyr)
library(stringr)

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate5_MPX <- read.csv("CKA_Multiplex_Fifth_24_qPCR_CqResults.csv")
dim (Plate5_MPX)

Plate5_Telo <- read.csv("CKA_Telomeres_Fifth_24_qPCR_CqResults.csv")
dim (Plate5_Telo)

# Concatenate data across runs for the same plate of samples
Plate5<- rbind(Plate5_MPX, Plate5_Telo)
dim(Plate5)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate5$PlateID <- "Plate5"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate5$Target[Plate5$Fluor == "Cy5"] <- "scnag"
Plate5$Target[Plate5$Fluor == "HEX"] <- "mtdna"
Plate5$Target[Plate5$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate5$Diff_AVG_Cq <- abs(Plate5$Cq - Plate5$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.3 threshold
Plate5$Flag_outlier <- ifelse(Plate5$Target %in% c("mtdna", "scnag") & Plate5$Diff_AVG_Cq > 0.3, "yes", "no")

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
  rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)

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

## ## Plate 6##
setwd("/Users/anetfilipova/Library/CloudStorage/Box-Box/TS_Lab_AnoleAging/NIH_R15_Anole_ModelAging_IGFs/Cox-Lab_Adult_KnownAges/AgeXSex_Telomeres/Data/Data_Analysis/Raw_Data_CqValues/Plate6")
library(dplyr)
library(stringr)

# Clear memory
rm(list=ls(all = TRUE))

# Import .csv files for each run for a particular plate of samples 
Plate6_MPX <- read.csv("CKA_Multiplex_Sixth_8_qPCR_CqResults.csv")
dim (Plate6_MPX)

Plate6_Telo <- read.csv("CKA_Telomeres_Sixth_8_qPCR_CqResults.csv")
dim (Plate6_Telo)

# Concatenate data across runs for the same plate of samples
Plate6<- rbind(Plate6_MPX, Plate6_Telo)
dim(Plate6)

# Add a column called PlateID and fill in correct Plate number to use as a variable in the statistics
Plate6$PlateID <- "Plate6"

# CHECK AND EDIT FOR YOUR DATA. 
#Name Correct Targets based on the fluorphores used in your reaction.
Plate6$Target[Plate6$Fluor == "Cy5"] <- "scnag"
Plate6$Target[Plate6$Fluor == "HEX"] <- "mtdna"
Plate6$Target[Plate6$Fluor == "SYBR"] <- "telomeres"

############## Look for outliers - 1st round ##############

# Calculate the absolute difference from the mean Cq
Plate6$Diff_AVG_Cq <- abs(Plate6$Cq - Plate6$Cq.Mean)

# Identify outliers for mtdna and scnag samples based on >0.3 threshold
Plate6$Flag_outlier <- ifelse(Plate6$Target %in% c("mtdna", "scnag") & Plate6$Diff_AVG_Cq > 0.3, "yes", "no")

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
  rename(Cq_Telomeres = Cq, Cq.Mean_Telomeres = Cq.Mean, Flag_outlier_Telomeres = Flag_outlier, SQ.Mean_Telomeres = SQ.Mean)

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
# Set working directory
setwd("/Users/anetfilipova/Library/CloudStorage/Box-Box/TS_Lab_AnoleAging/NIH_R15_Anole_ModelAging_IGFs/Cox-Lab_Adult_KnownAges/AgeXSex_Telomeres/Data/Data_Analysis/Raw_Data_CqValues")

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

#Merge
FinalData <- merge(qPCR_FinalData, Trait, by = c("Sample"))

# Save the merged data to a final CSV file
write.csv(FinalData, "CKA_FinalData.csv", row.names = FALSE)

# Optional: Print the dimensions of the final merged data
print(dim(FinalData))


## Running a linear regression model with Telomeres as dependent and Age_Months and Sex as independent variables
# Load data
data <- read.csv("CKA_FinalData.csv")
head(data)

# Install and load necessary packages
install.packages("lme4")
install.packages("lmerTest")
library(lme4)
library(lmerTest)

# Convert Plate_ID, Sex and Age_Class to factors and Age_Months to numeric
data$Plate_ID <- as.factor(data$Plate_ID)
data$Sex <- as.factor(data$Sex)
data$Age_Months <- as.numeric(data$Age_Months)
data$Age_Class <- as.factor(data$Age_Class)

######  Question 1. Do telomeres change with age, and is this different between the sexes?
telo_model <- lm(Telomeres.per.cell ~ Age_Months * Sex, data = data)

# Check the model summary
summary(telo_model)

data$Age_Class <- factor(data$Age_Months, levels = c("1", "7", "18", "36", "48", "60"))

# Summary table of remaining sample sizes using psych package
install.packages("psych")
library(psych)
summary_table <- table(FinalData$Age_Months, FinalData$Sex)
summary_table <- as.data.frame(summary_table)

# Install packages to visualize data
install.packages("emmeans")
install.packages("ggplot2")
install.packages("ggeffects")
install.packages("datawizard")
library(emmeans)
library(ggplot2)
library(ggeffects)
library(datawizard)


## Create a plot with trend lines and shaded confidence intervals to visualize the slopes for both sexes
p <- ggplot(FinalData, aes(x = Age_Months, y = Telomeres.per.cell, color = Sex, group = Sex)) +
  geom_point(size = 3, shape = 21, fill = "white") +  # Use shape 21 for filled points
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +  # Add linear trend lines with shaded confidence intervals
  labs(x = "Age (Months)", y = "Telomere repeats per cell") +  # Labels for axes
  scale_x_continuous(breaks = c(1, 7, 18, 36, 48, 60)) +  # Customize x-axis ticks
  theme_minimal()  # Optional: Customize the theme

## Print the plot
print(p)

## Save file as PNG for final figure production
ggsave(p, file="Telomere_Length_by_Age.png", width=9, height=7, dpi=600)

## Question 2. At any age, are telomeres different between sexes and between age classes?
# Fit a linear regression with Telomeres and Age_Class and Sex
telo_model_AgeClass <- lm(Telomeres.per.cell ~ Age_Class * Sex, data = data)

## Check the model summary
summary(telo_model_AgeClass)

emmeans(telo_model_AgeClass, pairwise ~ Age_Class * Sex, adjust = "tukey")
emmeans_data <- emmeans(telo_model_AgeClass, ~ Age_Class * Sex)

## Obtain estimated marginal means for Age_Months * Sex interaction
emm <- emmeans(telo_model_AgeClass, ~ Age_Class * Sex)

## Contrast for Sex differences within each Age_Class
contrast(emm, method = "pairwise", by = "Age_Class", adjust = "tukey")

## Contrast for Age_Class differences within each Sex
contrast(emm, method = "pairwise", by = "Sex", adjust = "tukey")

##Violin plot by sex
p2 = ggplot(FinalData, aes(x=Age_Class, y=Telomeres.per.cell, fill=Sex)) +
  ylab("Telomere Length") +
  xlab("Age (Months") +
  ggtitle("Telomere Length (T/S ratio)") +
  geom_violin(trim=F, position=position_dodge(width=0.9), scale="width", adjust = 0.60) +
  scale_fill_manual(values = c("red3", "skyblue3"), labels = c("Female", "Male")) + # Adjust colors as needed
  geom_boxplot(width = 0.3, color = "black", alpha = 0.5, position=position_dodge(width=0.9), outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.05, dodge.width = 0.9), size=1, alpha=0.5, color="black") +
  theme(axis.text.x = element_text(size=12), legend.position = "right", # Adjust legend position as needed
        plot.title = element_text(hjust = 0.5)) # Center the title
p2

# Save file as PNG for final figure production
ggsave(p2, file="Telomere_Length_by_Age.png", width=9, height=7, dpi=600)

######################
### Mitochondria - quick look
######  Question 1. Does mitochondrial DNA copy number change with age, and is this different between the sexes?
mtDNA_model <- lm(mtDNA.Mean ~ Age_Months * Sex, data = data)

# Check the model summary
summary(mtDNA_model)

# Get pairwise comparisons of slopes for Age_Months by Sex
emmeans(mtDNA_model, pairwise ~ Sex | Age_Months)

# Get estimated slopes for each sex
emmeans(mtDNA_model, specs = ~ Age_Months | Sex)

# Subset data for each sex
female_data <- subset(data, Sex == "F")
male_data <- subset(data, Sex == "M")

# Fit models for each sex separately
mtDNA_model_female <- lm(mtDNA.Mean ~ Age_Months, data = female_data)
mtDNA_model_male <- lm(mtDNA.Mean ~ Age_Months, data = male_data)

# Summarize each model
summary(mtDNA_model_female)
summary(mtDNA_model_male)


## Question 2. At any age, is the mtDNA copy number different between sexes and between age classes?
# Fit a linear regression with mtDNA as dependent and Age_Class and Sex as independent variables
mtDNA_model_Age_Class <- lm(mtDNA.Mean ~ Age_Class * Sex, data = data)
summary(mtDNA_model_Age_Class)

## Obtain estimated marginal means for Age_Class* Sex interaction
emm <- emmeans(mtDNA_model_Age_Class, ~ Age_Class * Sex)

## Contrast for Sex differences within each Age_Class
contrast(emm, method = "pairwise", by = "Age_Class", adjust = "tukey")

## Contrast for Age_Class differences within each Sex
contrast(emm, method = "pairwise", by = "Sex", adjust = "tukey")

## Create the plot with trend lines and shaded confidence intervals
p3 <- ggplot(FinalData, aes(x = Age_Months, y = mtDNA.Mean, color = Sex, group = Sex)) +
  geom_point(size = 3, shape = 21, fill = "white") +  # Use shape 21 for filled points
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +  # Add linear trend lines with shaded confidence intervals
  labs(x = "Age (Months)", y = "Mitochondrial DNA Copy Number per Blood Cell") +  # Labels for axes
  scale_x_continuous(breaks = c(1, 7, 18, 36, 48, 60)) +  # Customize x-axis ticks
  theme_minimal()  # Optional: Customize the theme

## Print the plot
print(p3)


# Save file as PNG for final figure production
ggsave(p3, file="Mitochondrial_DNA_Copy.png", width=9, height=7, dpi=600)

##Violin plot by mitochondrial by age and sex
p4 = ggplot(data, aes(x=Age_Class, y=mtDNA.Mean, fill=Sex)) +
  ylab("MtDNA") +
  xlab("Age (Months)") +
  ggtitle("MtDNA copy number per Blood Cell") +
  geom_violin(trim=F, position=position_dodge(width=0.9), scale="width", adjust = 0.60) +
  scale_fill_manual(values = c("red3", "skyblue3"), labels = c("Female", "Male")) + # Adjust colors as needed
  geom_boxplot(width = 0.3, color = "black", alpha = 0.5, position=position_dodge(width=0.9), outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.05, dodge.width = 0.9), size=1, alpha=0.5, color="black") +
  theme(axis.text.x = element_text(size=12), legend.position = "right", # Adjust legend position as needed
        plot.title = element_text(hjust = 0.5)) # Center the title
p4

# Save file as PNG for final figure production
ggsave(p4, file="MtDNA.png", width=9, height=7, dpi=600)



#Separate both sexes and test if the slopes for both sexes is indeed different, i.e. increasing for males and decreasing for females##

# Filter data for females
female_data <- subset(data, Sex == "F")

# Fit linear model for females
mtDNA_model_female <- lm(mtDNA.Mean ~ Age_Class, data = female_data)
summary(mtDNA_model_female)

# Filter data for males
male_data <- subset(data, Sex == "M")

# Fit linear model for males
mtDNA_model_male <- lm(mtDNA.Mean ~ Age_Class, data = male_data)
summary(mtDNA_model_male)

