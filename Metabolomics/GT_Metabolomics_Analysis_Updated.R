#Load relevant libraries

library(dplyr)
library(ggplot2)
library(lubridate)
library(lmodel2)
library(emmeans)
library(lme4)
library(data.table)
library(car)
library(Rmisc)
library(ggrepel)
library(gridExtra)


# Clear memory
rm(list=ls(all = TRUE))

setwd("C:/Users/aliam/Box/TS_Lab_GopherTortoises/Metabolomics/Data Analysis")

##### DATASET PRE_PROCESSING #####

#READ DATASETS & INSPECT STRUCTURE

Baseline_Bleed <- read.csv("C:/Users/aliam/Box/TS_Lab_GopherTortoises/Metabolomics/Data Analysis/Baseline_bleed_November_23.csv" ,header=T, sep = ",", as.is=T)
str(Baseline_Bleed)

Hibernation_Bleed <- read.csv("C:/Users/aliam/Box/TS_Lab_GopherTortoises/Metabolomics/Data Analysis/Hibernation_Bleed_2_27-28_24.csv" ,header=T, sep = ",", as.is=T)
str(Hibernation_Bleed)

Post_Hibernation_Bleed_2week <- read.csv("C:/Users/aliam/Box/TS_Lab_GopherTortoises/Metabolomics/Data Analysis/Post_hibernation_bleeding_03_19_24.csv" ,header=T, sep = ",", as.is=T)
str(Post_Hibernation_Bleed_2week)

Post_Hibernation_Bleed_3month <- read.csv("C:/Users/aliam/Box/TS_Lab_GopherTortoises/Metabolomics/Data Analysis/Post_Hibernation_3month_Bleed_06_09_24.csv" ,header=T, sep = ",", as.is=T)
str(Post_Hibernation_Bleed_3month)

AU_Headstart23_Individual_Data <- read.csv("C:/Users/aliam/Box/TS_Lab_GopherTortoises/Metabolomics/Data Analysis/Headstart_Individuals_AU_2023.csv" ,header=T, sep = ",", as.is=T)
str(AU_Headstart23_Individual_Data)

Glucose_Standard_Data <- read.csv("C:/Users/aliam/Box/TS_Lab_GopherTortoises/Metabolomics/Data Analysis/Glucose_Standard_Data.csv")
str(Glucose_Standard_Data)

TRIG_Master_Dataset <- read.csv("C:/Users/aliam/Box/TS_Lab_GopherTortoises/Metabolomics/Data Analysis/Triglyceride_MasterDataSheet_AU23GT.csv")
str(TRIG_Master_Dataset)

ECOA_Master_Dataset <- read.csv("C:/Users/aliam/Box/TS_Lab_GopherTortoises/Metabolomics/Data Analysis/ECOA_MasterDataSheet_AU23GT.csv")
str(ECOA_Master_Dataset)

# REMOVE UNNECESSARY COLUMNS

Updated_Baseline_Bleed <- Baseline_Bleed[, -c(4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 17)]
Updated_Baseline_Bleed

Updated_Hibernation_Bleed <- Hibernation_Bleed[, -c(4, 5, 6, 8, 9, 10, 11, 12, 13, 16)]
Updated_Hibernation_Bleed

Updated_Post_Hibernation_Bleed_2week <- Post_Hibernation_Bleed_2week[, -c(5, 6, 8, 9, 10, 11, 12, 13, 14, 17)]
Updated_Post_Hibernation_Bleed_2week

Updated_Post_Hibernation_Bleed_3month <- Post_Hibernation_Bleed_3month[, -c(5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 18)]
Updated_Post_Hibernation_Bleed_3month

Updated_AU_Headstart23_Individual_Data <- AU_Headstart23_Individual_Data[, -c(6,11)]
Updated_AU_Headstart23_Individual_Data

#Removing columns with "Not_Viable" Status (Individuals that failed to hatch)

Updated_AU_Headstart23_Individual_Data <- AU_Headstart23_Individual_Data %>%
  filter_all(all_vars(!grepl("Not_Viable", ., ignore.case = TRUE)))

Updated_AU_Headstart23_Individual_Data

standardize_types <- function(dfs) {
  all_colnames <- unique(unlist(lapply(dfs, colnames)))
  common_types <- sapply(all_colnames, function(col) {
    types <- unlist(lapply(dfs, function(df) if (col %in% colnames(df)) class(df[[col]]) else NA))
    types <- types[!is.na(types)]
    names(sort(table(types), decreasing = TRUE))[1]
  })
  lapply(dfs, function(df) {
    for (col in all_colnames) {
      df[[col]] <- if (col %in% colnames(df))
        switch(common_types[[col]],
               integer = as.integer(df[[col]]),
               numeric = as.numeric(df[[col]]),
               character = as.character(df[[col]]),
               factor = as.factor(df[[col]]),
               POSIXct = as.POSIXct(df[[col]]),
               df[[col]])
      else NA
    }
    df
  })
}
data_frames <- list(Updated_Baseline_Bleed,
                    Updated_Hibernation_Bleed,
                    Updated_Post_Hibernation_Bleed_2week,
                    Updated_Post_Hibernation_Bleed_3month
)
data_frames <- standardize_types(data_frames)
# Combine Data frames
Master_Bleed_Dataset <- bind_rows(data_frames)
Master_Bleed_Dataset <- Master_Bleed_Dataset %>%
  select(-which(duplicated(names(Master_Bleed_Dataset))))

str(Master_Bleed_Dataset)

# Add Treatment Column to Master_Dataset

Master_Bleed_Dataset <- Master_Bleed_Dataset %>%
  mutate(treatment = if_else(
    Clutch_ID %in% c(6.2, 9.2, 10.3, 9.4, 3.3, 8.1, 15.3, 10.1, 5.5, 11.4, 10.5, 2.4, 15.5, 12.6, 9.1, 9.5, 5.1, 15.1, 14.4, 11.2, 13.2, 6.1, 7.4, 7.2, 12.2, 5.3, 6.5, 17.3, 17.5, 13.4),
    "Constant Growth",
    "Cold Dormancy"
  ))

head(Master_Bleed_Dataset)

# Define key dates with mm/dd/yyyy format
baseline_dates <- as.Date(c("11/01/2023", "11/05/2023"), format = "%m/%d/%Y")
hibernation_dates <- as.Date(c("02/27/2024", "02/28/2024"), format = "%m/%d/%Y")
two_weeks_date <- as.Date("03/19/2024", format = "%m/%d/%Y")
three_months_date <- as.Date("06/10/2024", format = "%m/%d/%Y")

# Convert Date column to Date type with mm/dd/yyyy format
Master_Bleed_Dataset$Date <- as.Date(Master_Bleed_Dataset$Date, format = "%m/%d/%Y")


# Create a new column 'timepoint' based on the Date column
Master_Dataset <- Master_Bleed_Dataset %>%
  mutate(timepoint = case_when(
    Date %in% baseline_dates ~ "Pre-Dormancy",
    Date %in% hibernation_dates ~ "Dormancy",
    Date %in% two_weeks_date ~ "3 weeks-Post",
    Date %in% three_months_date ~ "3 months-Post",
    TRUE ~ "other"  
  ))

# Convert Date column to Date type with mm/dd/yyyy format
Master_Dataset$Date <- as.Date(Master_Dataset$Date, format = "%m/%d/%Y")

Master_Dataset$timepoint <- factor(Master_Dataset$timepoint, 
                                   levels = c("Pre-Dormancy", "Dormancy", "3 weeks-Post", "3 months-Post"))
print(Master_Dataset$timepoint)

### REMOVE INDIVIDUALS FROM EXPERIMENT

#Removed due to unidentified illness  
Master_Dataset <- Master_Dataset %>%
  filter(!Clutch_ID %in% c("10.3", "5.6", "6.1", "15.4"))

### Prepare TRIG dataset

TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  mutate(
    Timepoint = sub(".*\\((.*)\\)", "\\1", Clutch_ID), # Extract the letter in parentheses
    Clutch_ID = sub(" \\(.*\\)", "", Clutch_ID)       # Remove the parentheses and letter
  )

# Remove empty rows
TRIG_Master_Dataset <- na.omit(TRIG_Master_Dataset)
TRIG_Master_Dataset 

# Update TRIG values for Timepoint A where TRIG is negative
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  mutate(TRIG = if_else(Timepoint == "A" & TRIG < 0, 0, TRIG))

# Remove rows where TRIG is negative for all timepoints except A
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  filter(!(Timepoint != "A" & TRIG < 0))

#Remove rows where CV is > 10

TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  filter(CV <= 10)

#Remove rows where Plasma.Appearance is not U
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  filter(Plasma.Appearance == "U")

# Extract rows where Treatment is POOL into a new dataset
TRIG_Assay_Variation <- TRIG_Master_Dataset %>%
  filter(Treatment == "POOL")

# Remove rows where Treatment is POOL from TRIG_Master_Dataset
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  filter(Treatment != "POOL")

# Remove rows where Clutch_ID is NA
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  filter(Clutch_ID != "" & !is.na(Clutch_ID))

# Remove rows where Clutch_ID is blank
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  filter(Clutch_ID != "BLANK")

#Create log-transformed CORT and TRIG columns
TRIG_Master_Dataset$logTRIG <- log10(TRIG_Master_Dataset$TRIG)

# Replace -Inf values with 0 in the logTRIG column
TRIG_Master_Dataset$logTRIG[TRIG_Master_Dataset$logTRIG == -Inf] <- 0

# Identify (and print) Clutch_IDs that have a timepoint with exactly 2 TRIG entries
duplicate_clutch_ids <- TRIG_Master_Dataset %>%
  group_by(Clutch_ID, Timepoint) %>%
  filter(n() == 2) %>%         # keep only groups with exactly 2 rows
  ungroup() %>%
  distinct(Clutch_ID)          # extract unique Clutch_ID values

# Print the Clutch_IDs
print(duplicate_clutch_ids)


TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  # Add a count column for each Clutch_ID and Timepoint combination
  add_count(Clutch_ID, Timepoint, name = "dup_count") %>%
  # For groups with duplicates (dup_count > 1), filter out rows with the undesired Plate
  filter(!(dup_count > 1 & (
    (Clutch_ID == 12.5 & Plate == "3") |
      (Clutch_ID == 13.1 & Plate == "4") |
      (Clutch_ID == 13.4 & Plate == "4") |
      (Clutch_ID == 15.1 & Plate == "2") |
      (Clutch_ID == 17.5 & Plate == "3")
  ))) %>%
  # Remove the helper column
  select(-dup_count)

print(TRIG_Master_Dataset)

# Merge Tank & Nest_ID columns from Master Individual Data sheet to Triglycerides

TRIG_Master_Dataset <- merge(
  x = TRIG_Master_Dataset,
  y = Updated_AU_Headstart23_Individual_Data[, c("Clutch_ID", "Tank", "Nest_ID")],
  by = "Clutch_ID",
  all.x = FALSE  
)
str(TRIG_Master_Dataset)

# Remove individuals for bleed times > 10
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  filter(!Clutch_ID %in% c("9.3", "2.2", "14.4", "10.4", "11.2", "3.4"))

TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  mutate(Timepoint = trimws(Timepoint),         # Remove leading/trailing spaces
         Timepoint = toupper(Timepoint))       # Convert to uppercase for consistency

unique(TRIG_Master_Dataset$Timepoint)

# Calculate z-scores
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  mutate(z_score = (TRIG - mean(TRIG, na.rm = TRUE)) / sd(TRIG, na.rm = TRUE))

# Identify outliers with z-scores > 3 or < -3
outliers <- TRIG_Master_Dataset %>%
  filter(abs(z_score) > 3)

print(outliers)

# Filter and save data for Timepoint A
timepoint_A_data <- TRIG_Master_Dataset %>%
  filter(Timepoint == "A")

# Filter and save data for Timepoint C
timepoint_C_data <- TRIG_Master_Dataset %>%
  filter(Timepoint == "C")

# Filter and save data for Timepoint D
timepoint_D_data <- TRIG_Master_Dataset %>%
  filter(Timepoint == "D")

# Convert A, C, D in Timepoint column 
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  mutate(Timepoint = case_when(
    Timepoint == "A" ~ "Pre-Dormancy",
    Timepoint == "C" ~ "3 weeks (Post)",
    Timepoint == "D" ~ "3 months (Post)",
    TRUE ~ Timepoint  # Keeps other values unchanged
  ))

unique(TRIG_Master_Dataset$Timepoint)

# Filter individuals with values for all three timepoints
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  group_by(Clutch_ID) %>%
  filter(all(c("Pre-Dormancy", "3 weeks (Post)", "3 months (Post)") %in% Timepoint)) %>%
  ungroup()

# Reorder the levels of the Timepoint variable
TRIG_Master_Dataset$Timepoint <- factor(
  TRIG_Master_Dataset$Timepoint, 
  levels = c("Pre-Dormancy", "3 weeks (Post)", "3 months (Post)")
)

#Calculate within plate variation

within_plate_variation <- TRIG_Master_Dataset %>%
  group_by(Plate) %>%
  summarize(total_CV = mean(CV, na.rm = TRUE))

print(within_plate_variation)

## Definining an object for a colour-blind friendly pallette ##
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### TRIG ANALYSIS ####

# Calculate sample size per treatment
sample_size_per_treatment <- TRIG_Master_Dataset %>%
  group_by(Treatment) %>%
  summarise(Sample_Size = n_distinct(Clutch_ID))

# Print sample size per treatment
print(sample_size_per_treatment)


mean(TRIG_Master_Dataset$TRIG)
sd(TRIG_Master_Dataset$TRIG)
min(TRIG_Master_Dataset$TRIG)
max(TRIG_Master_Dataset$TRIG)
hist(TRIG_Master_Dataset$TRIG)
hist(TRIG_Master_Dataset$logTRIG)


TRIG_Model <- lme(TRIG ~  Treatment + Timepoint, 
                      random = ~1 | Nest_ID,
                      data = TRIG_Master_Dataset)

# Summary of the model
summary(TRIG_Model)
anova(TRIG_Model)
shapiro.test(residuals(TRIG_Model))
hist(residuals(TRIG_Model))

#TreatmentOverwinter:Timepoint3 weeks - Not significant (pvalue - 0.4339)
#TreatmentOverwinter:Timepoint3 months - Not significant (pvalue - 0.8210)


ggplot(TRIG_Master_Dataset, aes(x = Timepoint, y = TRIG, fill = Treatment)) +
  geom_boxplot() +
  labs(title = "TRIG Levels by Treatment & Timepoint",
       x = "Timepoint",
       y = "Triglyceride Concentration (mg/dL)") +
)  

pairwise_interaction <- emmeans(TRIG_Model, pairwise ~ Treatment * Timepoint, adjust = "tukey")
print(pairwise_interaction)

#Pre-Dormancy
TRIG_Model_Pre <- lme(TRIG ~  Treatment,
                     random = ~1 | Nest_ID,
                     data = timepoint_A_data)

# Summary of the model
summary(TRIG_Model_Pre)
anova(TRIG_Model_Pre)
shapiro.test(residuals(TRIG_Model_Pre))
hist(residuals(TRIG_Model_Pre))

plotA <- ggplot(timepoint_A_data, aes(x = Treatment, y = TRIG, color = Treatment)) +
  # Add boxplots with adjusted dodge to avoid overlap
  geom_boxplot(position = position_dodge(0.85)) +
  # Add jittered points for individual data visualization
  geom_jitter(width = 0.10, alpha = 0.5, size = 2) +
  # Label the axes
  ylab("Total Triglyceride Concentration (mg/dL)") +
  xlab("Treatment") +
  # Apply the same manual color scale (using specific indices from the palette)
  scale_color_manual(values = c(cbbPalette[6], cbbPalette[7]), 
                     name = "", 
                     labels = c("", "")) +
  # Use a classic theme and customize text sizes and remove legend
  theme_classic() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(color = "black", size = 14,face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 15))) +
  # Add the letter "C" at the top left (near the y-axis) as a figure label
  annotate("text", 
           x = -Inf, y = Inf, 
           label = "A", 
           hjust = -0.5, # Adjust horizontal justification to move it outside the plotting area
           vjust = 1.5,  # Adjust vertical justification to fine-tune its placement
           size = 8, 
           color = "black", 
           fontface = "bold")

plotA

#3week
TRIG_Model_3w <- lme(TRIG ~  Treatment,
                  random = ~1 | Nest_ID,
                  data = timepoint_C_data)

# Summary of the model
summary(TRIG_Model_3w)
shapiro.test(residuals(TRIG_Model_3w))
hist(residuals(TRIG_Model_3w))


# Calculate a y offset based on your data range
maxTRIG <- max(timepoint_C_data$TRIG, na.rm = TRUE)
y_offset <- diff(range(timepoint_C_data$TRIG, na.rm = TRUE)) * 0.05  # 5% above max
line_offset <- diff(range(timepoint_C_data$TRIG, na.rm = TRUE)) * 0.02  # 2% below the asterisk

plotC <- ggplot(timepoint_C_data, aes(x = Treatment, y = TRIG, color = Treatment)) +
  # Add boxplots with adjusted dodge to avoid overlap
  geom_boxplot(position = position_dodge(0.85)) +
  # Add jittered points for individual data visualization
  geom_jitter(width = 0.10, alpha = 0.5, size = 2) +
  # Label the axes
  ylab("Total Triglyceride Concentration (mg/dL)") +
  xlab("Treatment") +
  # Apply the same manual color scale (using specific indices from the palette)
  scale_color_manual(values = c(cbbPalette[6], cbbPalette[7]), 
                     name = "", 
                     labels = c("", "")) +
  # Use a classic theme and customize text sizes and remove legend
  theme_classic() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(color = "black", size = 14,face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 15))) +
  # Add an asterisk at the top right corner of the plot
  annotate("text", 
           label = "*", 
           x = 1.5, xend = 1.8, 
           y = maxTRIG + y_offset - line_offset, 
           yend = maxTRIG + y_offset - line_offset, 
           size = 12, 
           color = "black") +
  # Draw a horizontal segment (line) underneath the asterisk
  # Here, adjust x positions if you want the line wider or narrower.
  annotate("segment", 
           x = 1.2, xend = 1.8, 
           y = maxTRIG + y_offset - line_offset, 
           yend = maxTRIG + y_offset - line_offset, 
           color = "black", 
           size = 1) +
  # Add the letter "B" at the top left (near the y-axis) as a figure label
  annotate("text", 
           x = -Inf, y = Inf, 
           label = "B", 
           hjust = -0.5, # Adjust horizontal justification to move it outside the plotting area
           vjust = 1.5,  # Adjust vertical justification to fine-tune its placement
           size = 8, 
           color = "black", 
           fontface = "bold")

plotC

#3 month


TRIG_Model_3m <- lme(TRIG ~  Treatment,
                     random = ~1 | Nest_ID,
                     data = timepoint_D_data)

# Summary of the model
summary(TRIG_Model_3m)
shapiro.test(residuals(TRIG_Model_3m))
hist(residuals(TRIG_Model_3m))

plotD <- ggplot(timepoint_D_data, aes(x = Treatment, y = TRIG, color = Treatment)) +
  # Add boxplots with adjusted dodge to avoid overlap
  geom_boxplot(position = position_dodge(0.85)) +
  # Add jittered points for individual data visualization
  geom_jitter(width = 0.10, alpha = 0.5, size = 2) +
  # Label the axes
  ylab("Total Triglyceride Concentration (mg/dL)") +
  xlab("Treatment") +
  # Apply the same manual color scale (using specific indices from the palette)
  scale_color_manual(values = c(cbbPalette[6], cbbPalette[7]), 
                     name = "", 
                     labels = c("", "")) +
  # Use a classic theme and customize text sizes and remove legend
  theme_classic() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(color = "black", size = 14,face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 15))) +
  # Add the letter "C" at the top left (near the y-axis) as a figure label
  annotate("text", 
           x = -Inf, y = Inf, 
           label = "C", 
           hjust = -0.5, # Adjust horizontal justification to move it outside the plotting area
           vjust = 1.5,  # Adjust vertical justification to fine-tune its placement
           size = 8, 
           color = "black", 
           fontface = "bold")

plotD

Combined_TRIG_Plot <- grid.arrange(plotA, plotC, plotD, nrow = 1)

ggsave(Combined_TRIG_Plot, file="CombinedTRIG_Fig.jpg", dpi = 250, width = 15, height = 8)

#### ECOA ANALYSIS ####

ECOA_Master_Dataset <- ECOA_Master_Dataset %>%
  filter(CV <= 15)

#Removed due to negative ECOA value

ECOA_Master_Dataset <- ECOA_Master_Dataset %>%
  filter(ECOA >= 0)

# Remove empty rows & columns

ECOA_Master_Dataset <- ECOA_Master_Dataset %>%
  select(where(~ !all(is.na(.))))

ECOA_Master_Dataset <- na.omit(ECOA_Master_Dataset)
ECOA_Master_Dataset 

ECOA_Master_Dataset <- ECOA_Master_Dataset %>%
  mutate(
    Clutch_ID = sub("\\(.\\)", "", Clutch_ID),       # Remove the parentheses and letter       
    )

# Convert both Clutch_ID columns to character and remove extra white space
ECOA_Master_Dataset$Clutch_ID <- trimws(as.character(ECOA_Master_Dataset$Clutch_ID))
Updated_AU_Headstart23_Individual_Data$Clutch_ID <- trimws(as.character(Updated_AU_Headstart23_Individual_Data$Clutch_ID))

# Examine unique values in each data frame
print(unique(ECOA_Master_Dataset$Clutch_ID))
print(unique(Updated_AU_Headstart23_Individual_Data$Clutch_ID))

# Now perform the left join
ECOA_Master_Dataset <- merge(
  ECOA_Master_Dataset,
  Updated_AU_Headstart23_Individual_Data[, c("Clutch_ID", "Nest_ID")],
  by = "Clutch_ID"
)
str(ECOA_Master_Dataset)

# Remove individuals for bleed times > 10
ECOA_Master_Dataset <- ECOA_Master_Dataset %>%
  filter(!Clutch_ID %in% c("9.3", "2.2", "14.4", "10.4", "11.2", "3.4"))


ECOA_Model_3w <- lme(ECOA ~  Treatment,
                     random = ~1 | Nest_ID,
                     data = ECOA_Master_Dataset)

# Summary of the model
summary(ECOA_Model_3w)
anova(ECOA_Model_3w)
shapiro.test(residuals(ECOA_Model_3w))
hist(residuals(ECOA_Model_3w))



# Calculate a y offset based on your data range
maxECOA <- max(ECOA_Master_Dataset$ECOA, na.rm = TRUE)
y_offset_ECOA <- diff(range(ECOA_Master_Dataset$ECOA, na.rm = TRUE)) * 0.05  # 5% above max
line_offset_ECOA <- diff(range(ECOA_Master_Dataset$ECOA, na.rm = TRUE)) * 0.02  # 2% below the asterisk

ECOA_plot <- ggplot(ECOA_Master_Dataset, aes(x = Treatment, y = ECOA, color = Treatment)) +
  # Add boxplots with adjusted dodge to avoid overlap
  geom_boxplot(position = position_dodge(0.85)) +
  # Add jittered points for individual data visualization
  geom_jitter(width = 0.10, alpha = 0.5, size = 2) +
  # Label the axes
  ylab("Acetyl CoA Absorbances") +
  xlab("Timepoint") +
  # Apply the same manual color scale (using specific indices from the palette)
  scale_color_manual(values = c(cbbPalette[6], cbbPalette[7]), 
                     name = "", 
                     labels = c("", "")) +
  # Use a classic theme and customize text sizes and remove legend
  theme_classic() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(color = "black", size = 14,face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 15))) +
  # Add an asterisk at the top right corner of the plot
  annotate("text", 
           label = "*", 
           x = 1.5, xend = 1.8, 
           y = maxECOA + y_offset_ECOA - line_offset_ECOA, 
           yend = maxECOA + y_offset_ECOA - line_offset_ECOA, 
           size = 12, 
           color = "black") +
  # Draw a horizontal segment (line) underneath the asterisk
  # Here, adjust x positions if you want the line wider or narrower.
  annotate("segment", 
           x = 1.2, xend = 1.8, 
           y = maxECOA + y_offset_ECOA - line_offset_ECOA, 
           yend = maxECOA + y_offset_ECOA - line_offset_ECOA, 
           color = "black", 
           size = 1)

ECOA_plot

# Calculate sample size per treatment
sample_size_per_treatment_ECOA <- ECOA_Master_Dataset %>%
  group_by(Treatment) %>%
  summarise(Sample_Size = n_distinct(Clutch_ID))

# Print sample size per treatment
print(sample_size_per_treatment_ECOA)

within_plate_variation_ECOA <- ECOA_Master_Dataset %>%
  total_CV = mean(CV, na.rm = TRUE)

print(within_plate_variation_ECOA)

ggsave(ECOA_plot, file="ECOA_Fig.jpg", dpi = 250, width = 10, height = 8)


#### Glucose Analysis ####

# Create Glucose Dataset

Glucose_Data <- Master_Dataset
Glucose_Data

# Convert Glucose column to numeric
Glucose_Data$Glucose <- as.numeric(as.character(Glucose_Data$Glucose))

# Remove rows where Glucose is "-"
Glucose_Data <- Glucose_Data %>%
  filter(Glucose != "-")

# Filter individuals with values for all  timepoints
Glucose_Data <- Glucose_Data %>%
  group_by(Clutch_ID) %>%
  filter(all(c("Pre-Dormancy", "Dormancy", "3 weeks-Post", "3 months-Post") %in% timepoint)) %>%
  ungroup()

# Calculate sample size per treatment
sample_size_per_treatment <- Glucose_Data %>%
  group_by(treatment) %>%
  summarise(Sample_Size = n_distinct(Clutch_ID))

# Print sample size per treatment
print(sample_size_per_treatment)

# Reorder the levels of the Timepoint variable
Glucose_Data$timepoint <- factor(
  Glucose_Data$timepoint, 
  levels = c("Pre-Dormancy", "Dormancy", "3 weeks-Post", "3 months-Post")
)

# Reorder the levels of the treatment variable
Glucose_Data$treatment <- factor(
  Glucose_Data$treatment, 
  levels = c("Constant Growth", "Cold Dormancy")
)

Glucose_Model <- lme(Glucose ~  treatment * timepoint,
                  random = ~1 | Tank,
                  data = Glucose_Data)

# Summary of the model
summary(Glucose_Model)
anova(Glucose_Model)
shapiro.test(residuals(Glucose_Model))
hist(residuals(Glucose_Model))

Glucose_Plot <- ggplot(Glucose_Data, aes(x =timepoint, y = Glucose, fill = treatment)) +
  geom_boxplot() +
  labs(title = "Plasma Glucose Concentrations by Treatment & Timepoint",
       x = "Timepoint",
       y = "Plasma Glucose Concentrations") +
  common_layers

Glucose_Plot

ggsave(Glucose_Plot, file="Glucose_Fig.jpg", dpi = 250, width = 10, height = 8)

pairwise_interaction_gl <- emmeans(Glucose_Model, pairwise ~ treatment * timepoint, adjust = "tukey")
print(pairwise_interaction_gl)

# Filter and save data for Timepoint C
Glucose_timepoint_B_data <- Glucose_Data %>%
  filter(timepoint == "Dormancy")

# Filter and save data for Timepoint C
Glucose_timepoint_C_data <- Glucose_Data %>%
  filter(timepoint == "3 weeks-Post")

# Filter and save data for Timepoint D
Glucose_timepoint_D_data <- Glucose_Data %>%
  filter(timepoint == "3 months-Post")

Glucose_Model1 <- lme(Glucose ~  treatment,
                     random = ~1 | Clutch_ID,
                     data = Glucose_timepoint_B_data)

# Summary of the model
summary(Glucose_Model1)
anova(Glucose_Model1)
shapiro.test(residuals(Glucose_Model1))
hist(residuals(Glucose_Model1))

Glucose_Model2 <- lme(Glucose ~  treatment,
                      random = ~1 | Clutch_ID,
                      data = Glucose_timepoint_C_data)

# Summary of the model
summary(Glucose_Model2)
anova(Glucose_Model2)
shapiro.test(residuals(Glucose_Model2))
hist(residuals(Glucose_Model2))

Glucose_Model3 <- lme(Glucose ~  treatment,
                      random = ~1 | Clutch_ID,
                      data = Glucose_timepoint_D_data)

# Summary of the model
summary(Glucose_Model3)
anova(Glucose_Model3)
shapiro.test(residuals(Glucose_Model3))
hist(residuals(Glucose_Model3))

