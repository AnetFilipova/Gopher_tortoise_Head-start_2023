# Exploring the effects of first-year cold dormancy on metabolic profiles in head-started Gopher tortoises

#Load relevant libraries
library(dplyr)
library(ggplot2)
library(lubridate)
library(lmodel2)
library(emmeans)
library(lmerTest)
library(lme4)
library(data.table)
library(car)
library(Rmisc)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)


# Clear memory
rm(list=ls(all = TRUE))

##### DATASET PREPARATION #####

Baseline_Bleed <- read.csv("Metabolites/Data/Baseline_bleed_November_23.csv" ,header=T, sep = ",", as.is=T)
str(Baseline_Bleed)

Hibernation_Bleed <- read.csv("Metabolites/Data/Hibernation_Bleed_2_27-28_24.csv" ,header=T, sep = ",", as.is=T)
str(Hibernation_Bleed)

Post_Hibernation_Bleed_2week <- read.csv("Metabolites/Data/Post_hibernation_bleeding_03_19_24.csv" ,header=T, sep = ",", as.is=T)
str(Post_Hibernation_Bleed_2week)

Post_Hibernation_Bleed_3month <- read.csv("Metabolites/Data/Post_Hibernation_3month_Bleed_06_09_24.csv" ,header=T, sep = ",", as.is=T)
str(Post_Hibernation_Bleed_3month)

AU_Headstart23_Individual_Data <- read.csv("Metabolites/Data/Headstart_Individuals_AU_2023.csv" ,header=T, sep = ",", as.is=T)
str(AU_Headstart23_Individual_Data)

Glucose_Standard_Data <- read.csv("Metabolites/Data/Glucose_Standard_Data.csv")
str(Glucose_Standard_Data)

TRIG_Master_Dataset <- read.csv("Metabolites/Data/Triglyceride_MasterDataSheet_AU23GT.csv")
str(TRIG_Master_Dataset)

ECOA_Master_Dataset <- read.csv("Metabolites/Data/ECOA_MasterDataSheet_AU23GT.csv")
str(ECOA_Master_Dataset)

TRIG_Parallelism <- read.csv("Metabolites/Data/tg_dilution_linearity.csv")
str(TRIG_Parallelism)

# Demonstrate Parallelism of Triglyceride (do diluted standards match undiluted sample values)

## Definining an object for a colour-blind friendly pallette ##
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- ggplot( TRIG_Parallelism, aes(x = Concentration_mg_dL, y = Absorbance) ) + 
  # Triglyceride standards
  geom_point( 
    data = TRIG_Parallelism %>% filter(Type == "Standard"), 
    color = "black", 
    size = 3 
  ) + 
  geom_line( 
    data = TRIG_Parallelism %>% filter(Type == "Standard"), 
    color = "black", 
    linewidth = 1 
  ) + 
  # Plasma samples
  geom_point( 
    data = TRIG_Parallelism %>% filter(Type == "Sample"), 
    aes(color = Group), 
    size = 3 
  ) + 
  geom_line( 
    data = TRIG_Parallelism %>% filter(Type == "Sample"), 
    aes(color = Group, group = Group), 
    linewidth = 1 
  ) + 
  scale_color_manual( 
    values = c( 
      "Constant Warmth" = cbbPalette[7], 
      "Cold Dormancy" = cbbPalette[6] 
    ), 
    labels = c( 
      "Constant Warmth" = "Constant-Warmth", 
      "Cold Dormancy" = "Cold-Dormancy" 
    ),
    name = "Treatment" 
  ) +
  labs( 
    x = "Triglyceride Concentration (mg/dL)", 
    y = "Absorbance" 
  ) + 
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 20, b = 10)),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.title.align = 0.5
  )

print(p)

ggsave(filename = "Metabolites/Figures/Triglyceride_Parallelism.png", plot = p,width=9, height=7, dpi=600)

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
    "Constant-Warmth",
    "Cold-Dormancy"
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
  mutate(Timepoint = case_when(
    Date %in% baseline_dates ~ "Before dormancy",
    Date %in% hibernation_dates ~ "Dormancy end",
    Date %in% two_weeks_date ~ "3 weeks-Post",
    Date %in% three_months_date ~ "3 months-Post",
    TRUE ~ "other"  
  ))

Master_Dataset$Timepoint <- factor(Master_Dataset$Timepoint, 
                                   levels = c("Before dormancy", "Dormancy end", "3 weeks-Post", "3 months-Post"))
print(Master_Dataset$Timepoint)

### REMOVE INDIVIDUALS FROM EXPERIMENT

Master_Dataset <- Master_Dataset %>%
  filter(!Clutch_ID %in% c("5.6", "6.1", "15.4"))

############ TRIG ANALYSIS ###############

TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  mutate(
    Timepoint = sub(".*\\((.*)\\)", "\\1", Clutch_ID), # Extract the letter in parentheses
    Clutch_ID = sub(" \\(.*\\)", "", Clutch_ID)       # Remove the parentheses and letter
  )

# Remove empty rows
TRIG_Master_Dataset <- na.omit(TRIG_Master_Dataset)
TRIG_Master_Dataset 

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
      (Clutch_ID == 12.5 & Plate == "2") |
      (Clutch_ID == 9.3 & Plate == "3") |
      (Clutch_ID == 13.4 & Plate == "4")|
    (Clutch_ID == 13.1 & Plate == "2")|
      (Clutch_ID == 15.1 & Plate == "2") |
      (Clutch_ID == 17.5 & Plate == "3")
  ))) %>%
  # Remove the helper column
  select(-dup_count)

head(TRIG_Master_Dataset)

TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  group_by(Clutch_ID, Plate, Timepoint) %>%
  # within each (Clutch,Plate,Timepoint) group, drop the row
  # where CV is the maximum—but only when Clutch_ID==9.5 & Plate=="4"
  filter(!(Clutch_ID == 9.5 & Plate == "4" & CV == max(CV))) %>%
  ungroup()
print(TRIG_Master_Dataset)

# Check for leading/trailing spaces
TRIG_Master_Dataset$Clutch_ID <- trimws(TRIG_Master_Dataset$Clutch_ID)
TRIG_Master_Dataset$Timepoint <- trimws(TRIG_Master_Dataset$Timepoint)

Updated_AU_Headstart23_Individual_Data$Clutch_ID <- trimws(Updated_AU_Headstart23_Individual_Data$Clutch_ID)
Master_Dataset$Clutch_ID <- trimws(Master_Dataset$Clutch_ID)
Master_Dataset$Timepoint <- trimws(Master_Dataset$Timepoint)

# Merge Tank & Nest_ID columns from Master Individual Data sheet to Triglycerides

TRIG_Master_Dataset <- merge(
  x = TRIG_Master_Dataset,
  y = Updated_AU_Headstart23_Individual_Data[, c("Clutch_ID", "Tank", "Nest_ID")],
  by = "Clutch_ID",
  all.x = FALSE  
)
str(TRIG_Master_Dataset)

# Convert A, C, D in Timepoint column 
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  mutate(Timepoint = case_when(
    Timepoint == "A" ~ "Before dormancy",
    Timepoint == "C" ~ "3 weeks-Post",
    Timepoint == "D" ~ "3 months-Post",
    TRUE ~ Timepoint  # Keeps other values unchanged
  ))

str(TRIG_Master_Dataset)

#merge bleed times 
TRIG_Master_Dataset <- merge(
  TRIG_Master_Dataset,
  unique(Master_Dataset[, c("Clutch_ID", "Timepoint", "Total_bleed_time")]),
  by = c("Clutch_ID", "Timepoint")
)

#Convert bleed time column to MM:SS format
TRIG_Master_Dataset <- TRIG_Master_Dataset %>%
  mutate(
    # Standardize format: add ":00" if no colon present
    Total_bleed_time_clean = ifelse(
      grepl(":", Total_bleed_time),
      Total_bleed_time,
      paste0(Total_bleed_time, ":00")
    ),
    # Convert to hms format (assumes MM:SS)
    Total_bleed_time = ms(Total_bleed_time_clean)
  ) %>%
  select(-Total_bleed_time_clean) %>%  # Remove temporary column
  filter(Total_bleed_time <= ms("10:00"))  # 10 minutes

# Filter and save data for Timepoint A
timepoint_A_data <- TRIG_Master_Dataset %>%
  filter(Timepoint == "Before")

# Filter and save data for Timepoint C
timepoint_C_data <- TRIG_Master_Dataset %>%
  filter(Timepoint == "3-Week Post-Dormancy")

# Filter and save data for Timepoint D
timepoint_D_data <- TRIG_Master_Dataset %>%
  filter(Timepoint == "3-Month Post-Dormancy")

#Calculate within plate variation

within_plate_variation <- TRIG_Master_Dataset %>%
  group_by(Plate) %>%
  summarize(total_CV = mean(CV, na.rm = TRUE))

print(within_plate_variation)

#Convert TRIG to numeric 

TRIG_Master_Dataset$TRIG <- as.numeric(TRIG_Master_Dataset$TRIG)
str(TRIG_Master_Dataset)

## Definining an object for a colour-blind friendly pallette ##
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### TRIG ANALYSIS - Repeated measures ####

# Check the structure of the data
str(TRIG_Master_Dataset)

# Summary statistics
mean(TRIG_Master_Dataset$TRIG)
sd(TRIG_Master_Dataset$TRIG)
min(TRIG_Master_Dataset$TRIG)
max(TRIG_Master_Dataset$TRIG)
hist(TRIG_Master_Dataset$TRIG)

# Check sample sizes per Treatment and Timepoint combination
sample_size_summary <- TRIG_Master_Dataset %>%
  group_by(Treatment, Timepoint) %>%
  tally(name = "Sample_Size") %>%  
  ungroup()

print(sample_size_summary)

#### Fit Mixed Model ####
TRIG_model <- lmer(TRIG ~ Treatment * Timepoint + (1 | Nest_ID) + (1 | Clutch_ID), 
                   data = TRIG_Master_Dataset,
                   control = lmerControl(optimizer = "bobyqa"))

# Model summary
summary(TRIG_model)

# Type III ANOVA
anova(TRIG_model, type = 3)

# Get confidence intervals
confint(TRIG_model)

##############################################################################
#### Model Diagnostics ####
##############################################################################

# 1. NORMALITY OF RESIDUALS
par(mfrow = c(1, 2))
qqnorm(resid(TRIG_model), main = "Q-Q Plot: Residuals")
qqline(resid(TRIG_model), col = "red")
hist(resid(TRIG_model), main = "Histogram of Residuals", xlab = "Residuals", breaks = 20)
par(mfrow = c(1, 1))
shapiro.test(resid(TRIG_model))

# 2. HOMOGENEITY OF VARIANCE
plot(fitted(TRIG_model), resid(TRIG_model),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lwd = 2)

# Plot residuals by treatment (using boxplot for factor variables)
boxplot(resid(TRIG_model) ~ TRIG_Master_Dataset$Treatment,
        xlab = "Treatment", ylab = "Residuals",
        main = "Residuals by Treatment")
abline(h = 0, col = "red", lwd = 2, lty = 2)

leveneTest(resid(TRIG_model) ~ Treatment, data = TRIG_Master_Dataset)

# 3. OUTLIERS
std_resid <- scale(resid(TRIG_model))
plot(std_resid, main = "Standardized Residuals", 
     ylab = "Standardized Residuals", xlab = "Observation")
abline(h = c(-3, 3), col = "red", lty = 2)

outliers <- which(abs(std_resid) > 3)
if(length(outliers) > 0) {
  print(paste("Potential outliers at observations:", paste(outliers, collapse = ", ")))
  print(TRIG_Master_Dataset[outliers, c("Clutch_ID", "Treatment", "Timepoint", "TRIG")])
} else {
  print("No outliers detected (|standardized residual| > 3)")
}

##############################################################################
#### Post-hoc Comparisons ####
##############################################################################

emm_timepoint <- emmeans(TRIG_model, ~ Timepoint)
print(emm_timepoint)

# Pairwise comparisons between timepoints
pairs_timepoint <- pairs(emm_timepoint)
print(pairs_timepoint)

# Timepoint differences within each Treatment
emm_timepoint_by_treatment <- emmeans(TRIG_model, ~ Timepoint | Treatment)
print("Estimated marginal means by Timepoint within each Treatment:")
print(emm_timepoint_by_treatment)

# Pairwise comparisons 
pairs_cold <- pairs(emm_timepoint_by_treatment)
print("Pairwise timepoint comparisons:")
print(pairs_cold)

# Treatment differences at each timepoint
emm_treatment_by_timepoint <- emmeans(TRIG_model, ~ Treatment | Timepoint)
print("Estimated marginal means by Treatment within each Timepoint:")
print(emm_treatment_by_timepoint)

pairs_treatment <- pairs(emm_treatment_by_timepoint)
print("Pairwise treatment comparisons at each timepoint:")
print(pairs_treatment)

##############################################################################
#### Visualizations ####
##############################################################################

## Defining an object for a colour-blind friendly palette ##
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Prepare the TRIG data with proper timepoint labels
TRIG_Plot_Data <- TRIG_Master_Dataset %>%
  mutate(
    Timepoint_Label = factor(
      Timepoint,
      levels = c("Before dormancy", "3 weeks-Post", "3 months-Post"),
      labels = c("BD", "3wkPD", "3moPD")
    )
  )

# Separate Before (Pre-Dormancy) from other timepoints
pre_dormancy <- TRIG_Plot_Data %>% filter(Timepoint_Label == "BD")
post_dormancy <- TRIG_Plot_Data %>% filter(Timepoint_Label != "BD")  

# Define y positions for significance brackets
ypos <- max(TRIG_Plot_Data$TRIG, na.rm = TRUE) + 5

TRIG_ALL_Before_Combined_Plot <- ggplot() +
  # Pre-Dormancy - single gray box (combined treatments)
  geom_boxplot(data = pre_dormancy, aes(x = Timepoint_Label, y = TRIG),
               color = "black", fill = "gray70", outlier.shape = NA, 
               width = 0.6, linewidth = 0.6) +
  geom_point(data = pre_dormancy, aes(x = Timepoint_Label, y = TRIG),
             position = position_jitter(width = 0.1), 
             size = 2, color = "black", fill = "white", shape = 21, stroke = 0.5) +
  
  # Post-Dormancy timepoints - separate boxes by Treatment (white fill, colored outlines)
  geom_boxplot(data = post_dormancy, aes(x = Timepoint_Label, y = TRIG, color = Treatment),
               fill = "white", outlier.shape = NA, width = 0.6, 
               position = position_dodge(0.8), linewidth = 0.6) +
  geom_jitter(data = post_dormancy, aes(x = Timepoint_Label, y = TRIG, color = Treatment),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
              size = 2, alpha = 0.4) +
  
  # Mean points for post-dormancy timepoints only
  stat_summary(data = post_dormancy, 
               aes(x = Timepoint_Label, y = TRIG, color = Treatment, group = Treatment),
               fun = mean, geom = "point", size = 2,
               position = position_dodge(width = 0.8)) +
  
  # Mean point for Before dormancy (black, single combined point)
  stat_summary(data = pre_dormancy, 
               aes(x = Timepoint_Label, y = TRIG),
               fun = mean, geom = "point", size = 2, color = "black") +
  
  # Significance bracket: Cold-Dormancy (3 weeks to 3 months)
  # Left vertical line
  annotate("segment", x = 1.8, xend = 1.8, y = ypos + 10, yend = ypos + 20,
           color = cbbPalette[6], linewidth = 0.8) +
  # Horizontal line
  annotate("segment", x = 1.8, xend = 2.8, y = ypos + 20, yend = ypos + 20,
           color = cbbPalette[6], linewidth = 0.8) +
  # Right vertical line
  annotate("segment", x = 2.8, xend = 2.8, y = ypos + 10, yend = ypos + 20,
           color = cbbPalette[6], linewidth = 0.8) +
  # Asterisk
  annotate("text", x = 2.3, y = ypos + 25, label = "*", 
           size = 8, color = cbbPalette[6]) +
  
  # Significance bracket: Constant-Warmth (3 weeks to 3 months) - slightly higher
  # Left vertical line
  annotate("segment", x = 2.2, xend = 2.2, y = ypos + 35, yend = ypos + 45,
           color = cbbPalette[7], linewidth = 0.8) +
  # Horizontal line
  annotate("segment", x = 2.2, xend = 3.2, y = ypos + 45, yend = ypos + 45,
           color = cbbPalette[7], linewidth = 0.8) +
  # Right vertical line
  annotate("segment", x = 3.2, xend = 3.2, y = ypos + 35, yend = ypos + 45,
           color = cbbPalette[7], linewidth = 0.8) +
  # Asterisk
  annotate("text", x = 2.7, y = ypos + 50, label = "*", 
           size = 8, color = cbbPalette[7]) +
  
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
  
  # Colors matching your colorblind palette
  scale_color_manual(
    values = c("Cold-Dormancy" = cbbPalette[6], "Constant-Warmth" = cbbPalette[7]),
    name = "Treatment",
    labels = c("Cold-Dormancy", "Constant-Warmth")
  ) +
  
  labs(
    x = "Timepoint",
    y = "Total Triglyceride Concentration (mg/dL)"
  )

TRIG_ALL_Before_Combined_Plot

ggsave(TRIG_ALL_Before_Combined_Plot, file="Metabolites/Figures/Triglycerides.png", width=9, height=7, dpi=600)


####################### ECOA ANALYSIS #####################

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
  Updated_AU_Headstart23_Individual_Data[, c("Clutch_ID", "Nest_ID", "Tank")],
  by = "Clutch_ID"
)
str(ECOA_Master_Dataset)

# Remove individuals for bleed times > 10
ECOA_Master_Dataset <- ECOA_Master_Dataset %>%
  filter(!Clutch_ID %in% c("9.3", "2.2", "14.4", "10.4", "11.2", "3.4"))


#### Fit Mixed Model ####
ECOA_Model_3w <- lmer(ECOA ~ Treatment + (1|Tank), data = ECOA_Master_Dataset, control = lmerControl(optimizer = "bobyqa"))

# Model summary
summary(ECOA_Model_3w)

# Type III ANOVA
anova(ECOA_Model_3w, type = 3)

# Get confidence intervals
confint(ECOA_Model_3w)

##############################################################################
#### Model Diagnostics ####
##############################################################################

# 1. NORMALITY OF RESIDUALS
par(mfrow = c(1, 2))
qqnorm(resid(ECOA_Model_3w), main = "Q-Q Plot: Residuals")
qqline(resid(ECOA_Model_3w), col = "red")
hist(resid(ECOA_Model_3w), main = "Histogram of Residuals", xlab = "Residuals", breaks = 20)
par(mfrow = c(1, 1))
shapiro.test(resid(ECOA_Model_3w))

# 2. HOMOGENEITY OF VARIANCE
resid_df <- data.frame(
  residuals = resid(ECOA_Model_3w),
  fitted = fitted(ECOA_Model_3w),
  treatment = ECOA_Master_Dataset$Treatment[as.numeric(names(resid(ECOA_Model_3w)))]
)


# 2. HOMOGENEITY OF VARIANCE
plot(fitted(ECOA_Model_3w), resid(ECOA_Model_3w),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lwd = 2)

leveneTest(resid(ECOA_Model_3w) ~ Treatment, data = ECOA_Master_Dataset)

# 3. OUTLIERS
std_resid <- scale(resid(ECOA_Model_3w))
plot(std_resid, main = "Standardized Residuals", 
     ylab = "Standardized Residuals", xlab = "Observation")
abline(h = c(-3, 3), col = "red", lty = 2)

outliers <- which(abs(std_resid) > 3)
if(length(outliers) > 0) {
  print(paste("Potential outliers at observations:", paste(outliers, collapse = ", ")))
  print(ECOA_Master_Dataset[outliers, c("Clutch_ID", "Treatment", "ECOA")])
} else {
  print("No outliers detected (|standardized residual| > 3)")
}

##############################################################################
#### Post-hoc Comparisons ####
##############################################################################

# Treatment effect
emm_treatment <- emmeans(ECOA_Model_3w, ~ Treatment)
print(emm_treatment)

# Pairwise comparison
pairs_treatment <- pairs(emm_treatment)
print(pairs_treatment)

confint(ECOA_Model_3w)


### The 95 % C.I. for Cold-Dormancy is ± 0.17 (0.231-0.0565)
### The 95 % C.I. for Constant-Warmth is ± 0.2 (0.430-0.2274)

##############################################################################
#### Visualizations ####
##############################################################################

## Colour-blind friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Calculate y position for significance line
ypos <- max(ECOA_Master_Dataset$ECOA, na.rm = TRUE) + 0.02

# Boxplot comparing treatments
ECOA_Plot <- ggplot(ECOA_Master_Dataset, aes(x = Treatment, y = ECOA, color = Treatment)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    fill = "white",
    alpha = 0.5,
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.1,
    alpha = 0.4,
    size = 2
  ) +
  stat_summary(
    fun = mean, 
    geom = "point", 
    size = 3,
    shape = 16  # Filled circle
  ) +
  # Horizontal line only
  annotate("segment", x = 1, xend = 2, y = ypos, yend = ypos,
           color = "black", size = 0.8) +
  # Asterisk
  annotate("text", x = 1.5, y = ypos + 0.01, label = "*", 
           size = 8, color = "black") +
  
  scale_color_manual(
    values = c("Cold-Dormancy" = cbbPalette[6], "Constant-Warmth" = cbbPalette[7]),
    name = "Treatment",
    labels = c("Cold-Dormancy", "Constant-Warmth")
  ) +
  labs(
    x = "Treatment group",
    y = "Acetyl-CoA (AU)"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20, b = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "none",
    legend.title.align = 0.5
  )

ECOA_Plot

ggsave(ECOA_Plot, file="Metabolites/Figures/Acetyl-CoA.png", width=9, height=7, dpi=600)

# Calculate sample size per treatment
sample_size_per_treatment_ECOA <- ECOA_Master_Dataset %>%
  group_by(Treatment) %>%
  summarise(Sample_Size = n_distinct(Clutch_ID))

# Print sample size per treatment
print(sample_size_per_treatment_ECOA)

within_plate_variation_ECOA <- ECOA_Master_Dataset %>%
  summarise(total_CV = mean(CV, na.rm = TRUE))


print(within_plate_variation_ECOA)


####################### GLUCOSE Analysis #########################

# Create Glucose Dataset

Glucose_Data <- Master_Dataset
Glucose_Data


#Convert bleed time column to MM:SS format
Glucose_Data <- Glucose_Data %>%
  mutate(
    # Standardize format: add ":00" if no colon present
    Total_bleed_time_clean = ifelse(
      grepl(":", Total_bleed_time),
      Total_bleed_time,
      paste0(Total_bleed_time, ":00")
    ),
    # Convert to hms format (assumes MM:SS)
    Total_bleed_time = ms(Total_bleed_time_clean)
  ) %>%
  select(-Total_bleed_time_clean) %>%  # Remove temporary column
  filter(Total_bleed_time <= ms("10:00"))  # 10 minutes



# Convert Glucose column to numeric
Glucose_Data$Glucose <- as.numeric(as.character(Glucose_Data$Glucose))

# Remove rows where Glucose is "-"
Glucose_Data <- Glucose_Data %>%
  filter(Glucose != "-")

# Calculate sample size per treatment
sample_size_per_treatment <- Glucose_Data %>%
  group_by(treatment) %>%
  summarise(Sample_Size = n_distinct(Clutch_ID))

# Print sample size per treatment
print(sample_size_per_treatment)

# Reorder the levels of the Timepoint variable
Glucose_Data$Timepoint <- factor(
  Glucose_Data$Timepoint, 
  levels = c("Before dormancy", "Dormancy end", "3 weeks-Post", "3 months-Post")
)

str(Glucose_Data)

# Reorder the levels of the treatment variable
Glucose_Data$treatment <- factor(
  Glucose_Data$treatment, 
  levels = c("Cold-Dormancy", "Constant-Warmth")
)

# Now perform the left join
Glucose_Data <- merge(
  Glucose_Data,
  Updated_AU_Headstart23_Individual_Data[, c("Clutch_ID", "Nest_ID")],
  by = "Clutch_ID"
)

# Check the structure of the data
str(Glucose_Data)

# Summary statistics
mean(Glucose_Data$Glucose)
sd(Glucose_Data$Glucose)
min(Glucose_Data$Glucose)
max(Glucose_Data$Glucose)
hist(Glucose_Data$Glucose)

# Check sample sizes per Treatment and Timepoint combination
sample_size_summary <- Glucose_Data %>%
  group_by(treatment, Timepoint) %>%
  tally(name = "Sample_Size") %>%  
  ungroup()

print(sample_size_summary)

#### Fit Mixed Model ####
Glucose_model <- lmer(Glucose ~ treatment * Timepoint + (1 | Nest_ID) +(1 | Tank) + (1 | Clutch_ID), 
                   data = Glucose_Data,
                   control = lmerControl(optimizer = "bobyqa"))

# Model summary
summary(Glucose_model)

# Type III ANOVA
anova(Glucose_model, type = 3)

##############################################################################
#### Model Diagnostics ####
##############################################################################

# 1. NORMALITY OF RESIDUALS
par(mfrow = c(1, 2))
qqnorm(resid(Glucose_model), main = "Q-Q Plot: Residuals")
qqline(resid(Glucose_model), col = "red")
hist(resid(Glucose_model), main = "Histogram of Residuals", xlab = "Residuals", breaks = 20)
par(mfrow = c(1, 1))
shapiro.test(resid(Glucose_model))

# 2. HOMOGENEITY OF VARIANCE
plot(fitted(Glucose_model), resid(Glucose_model),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lwd = 2)

# Plot residuals by treatment (using boxplot for factor variables)
boxplot(resid(Glucose_model) ~ Glucose_Data$treatment,
        xlab = "Treatment", ylab = "Residuals",
        main = "Residuals by Treatment")
abline(h = 0, col = "red", lwd = 2, lty = 2)

leveneTest(resid(Glucose_model) ~ treatment, data = Glucose_Data)

# 3. OUTLIERS
std_resid <- scale(resid(Glucose_model))
plot(std_resid, main = "Standardized Residuals", 
     ylab = "Standardized Residuals", xlab = "Observation")
abline(h = c(-3, 3), col = "red", lty = 2)

outliers <- which(abs(std_resid) > 3)
if(length(outliers) > 0) {
  print(paste("Potential outliers at observations:", paste(outliers, collapse = ", ")))
  print(Glucose_Data[outliers, c("Clutch_ID", "treatment", "Timepoint", "Glucose")])
} else {
  print("No outliers detected (|standardized residual| > 3)")
}

######## Post-hoc Comparisons ########

# Estimated marginal means by timepoint
emm_timepoint <- emmeans(Glucose_model, ~ Timepoint)
print(emm_timepoint)

# Pairwise comparisons between timepoints
pairs_timepoint <- pairs(emm_timepoint)
print(pairs_timepoint)

# Timepoint differences within each Treatment
emm_timepoint_by_treatment <- emmeans(Glucose_model, ~ Timepoint | treatment)
print("Estimated marginal means by Timepoint within each Treatment:")
print(emm_timepoint_by_treatment)

# Pairwise comparisons within each treatment
pairs_by_treatment <- pairs(emm_timepoint_by_treatment)
print("Pairwise timepoint comparisons:")
print(pairs_by_treatment)

# Treatment differences within each Timepoint
emm_treatment_by_timepoint <- emmeans(Glucose_model, ~ treatment | Timepoint)
print("Estimated marginal means by Treatment within each Timepoint:")
print(emm_treatment_by_timepoint)

# Pairwise treatment comparisons at each timepoint
pairs_treatment_by_timepoint <- pairs(emm_treatment_by_timepoint)
print("Treatment comparisons at each timepoint:")
print(pairs_treatment_by_timepoint)

### The 95 % C.I. for End-Dormancy in Cold-Dormancy is ± 21.8 (52.8-31.0)
### The 95 % C.I. for End-Dormancy in Constant-Warmth is ± 21.6 (117.7-96.1)

############ VISUALIZATIONS ################


# Prepare the Glucose data with proper timepoint labels
Glucose_Plot_Data <- Glucose_Data %>%
  mutate(
    Timepoint_Label = factor(
      Timepoint,
      levels = c("Before dormancy", "Dormancy end", "3 weeks-Post", "3 months-Post"),
      labels = c("BD", "ED", "3wkPD", "3moPD")  
    )
  )

# Separate Before (Pre-Dormancy) from other timepoints
# IMPORTANT: Filter using the NEW label "Before", not the old level
pre_dormancy <- Glucose_Plot_Data %>% filter(Timepoint_Label == "BD")
post_dormancy <- Glucose_Plot_Data %>% filter(Timepoint_Label != "BD")

# Define y positions for significance brackets
ypos <- max(Glucose_Plot_Data$Glucose, na.rm = TRUE) + 10

Glucose_Plot <- ggplot() +
  # Pre-Dormancy - single gray box (combined treatments)
  geom_boxplot(data = pre_dormancy, aes(x = Timepoint_Label, y = Glucose),
               color = "black", fill = "gray70", outlier.shape = NA, 
               width = 0.6, linewidth = 0.6) +
  geom_point(data = pre_dormancy, aes(x = Timepoint_Label, y = Glucose),
             position = position_jitter(width = 0.1), 
             size = 2, color = "black", fill = "white", shape = 21, stroke = 0.5) +
  
  # Post-Dormancy timepoints - separate boxes by Treatment (white fill, colored outlines)
  geom_boxplot(data = post_dormancy, aes(x = Timepoint_Label, y = Glucose, color = treatment),
               fill = "white", outlier.shape = NA, width = 0.6, 
               position = position_dodge(0.8), linewidth = 0.6) +
  geom_jitter(data = post_dormancy, aes(x = Timepoint_Label, y = Glucose, color = treatment),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
              size = 2, alpha = 0.4) +
  
  # Mean points for post-dormancy timepoints only
  stat_summary(data = post_dormancy, 
               aes(x = Timepoint_Label, y = Glucose, color = treatment, group = treatment),
               fun = mean, geom = "point", size = 2,
               position = position_dodge(width = 0.8)) +
  
  # Mean point for Before dormancy (black, single combined point)
  stat_summary(data = pre_dormancy, 
               aes(x = Timepoint_Label, y = Glucose),
               fun = mean, geom = "point", size = 2, color = "black") +
  # Positioned just below the blue brackets
  annotate("segment", x = 1.8, xend = 2.2, y = ypos - 10, yend = ypos - 10,
           color = "black", linewidth = 0.8) +
  annotate("text", x = 2.0, y = ypos - 5, label = "****", 
           size = 6, color = "black") +
  
  # BLUE BRACKET 1: Cold-Dormancy Dormancy end → 3 weeks-Post (p < 0.001)
  annotate("segment", x = 1.8, xend = 1.8, y = ypos, yend = ypos + 10,
           color = cbbPalette[6], linewidth = 0.8) +
  annotate("segment", x = 1.8, xend = 2.8, y = ypos + 10, yend = ypos + 10,
           color = cbbPalette[6], linewidth = 0.8) +
  annotate("segment", x = 2.8, xend = 2.8, y = ypos, yend = ypos + 10,
           color = cbbPalette[6], linewidth = 0.8) +
  annotate("text", x = 2.3, y = ypos + 15, label = "****", 
           size = 6, color = cbbPalette[6]) +
  
  # BLUE BRACKET 2: Cold-Dormancy Dormancy end → 3 months-Post (p < 0.001)
  annotate("segment", x = 1.8, xend = 1.8, y = ypos + 20, yend = ypos + 30,
           color = cbbPalette[6], linewidth = 0.8) +
  annotate("segment", x = 1.8, xend = 3.8, y = ypos + 30, yend = ypos + 30,
           color = cbbPalette[6], linewidth = 0.8) +
  annotate("segment", x = 3.8, xend = 3.8, y = ypos + 20, yend = ypos + 30,
           color = cbbPalette[6], linewidth = 0.8) +
  annotate("text", x = 2.8, y = ypos + 35, label = "****", 
           size = 6, color = cbbPalette[6]) +
  
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
  
  # Colors matching your colorblind palette
  scale_color_manual(
    values = c("Cold-Dormancy" = cbbPalette[6], "Constant-Warmth" = cbbPalette[7]),
    name = "Treatment",
    labels = c("Cold-Dormancy", "Constant-Warmth")
  ) +
  
  labs(
    x = "Timepoint",
    y = "Plasma Glucose Concentration (mg/dL)"
  )

Glucose_Plot

ggsave(Glucose_Plot, file="Metabolites/Figures/Glucose.png", width=9, height=7, dpi=600)



####### Combining the 4 separate figures for mtDNA_Treatment, Glucose, Triglycerides and Acetyl-CoA into one 4-panel figure
## this object comes from the Telo_mtDNA script so we're importing it
mtdna_plot <- readRDS("Telo_mtDNA/Figures/mtdna_plot.rds")

# Extract single legend with transparent background
single_legend <- get_legend(
  mtdna_plot + 
    theme(legend.position = "right",
          legend.title = element_text(size = 17, face = "bold"),
          legend.text = element_text(size = 16),
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_blank())
)

# Remove all legends from individual plots and add tags outside panels
p1 <- mtdna_plot +
  labs(tag = "A") +
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold", size = 20),
        plot.tag.position = c(0, 1),
        plot.margin = margin(5, 30, 5, 5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(margin = margin(r = 3)))

p2 <- Glucose_Plot +
  labs(tag = "B") +
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold", size = 20),
        plot.tag.position = c(-0.04, 1),
        plot.margin = margin(5, 5, 5, 30),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(margin = margin(r = 3)))

p3 <- TRIG_ALL_Before_Combined_Plot +
  labs(tag = "C") +
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold", size = 20),
        plot.tag.position = c(0, 1.05),
        plot.margin = margin(5, 30, 5, 5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(margin = margin(r = 3)))

p4 <- ECOA_Plot +
  labs(tag = "D") +
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold", size = 20),
        plot.tag.position = c(-0.04, 1.05),
        plot.margin = margin(5, 5, 5, 30),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(margin = margin(r = 3)))

grid_no_legend <- (p1 + p2) / (p3 + p4)

combined_4panel <- plot_grid(grid_no_legend, single_legend,
                             ncol = 2,
                             rel_widths = c(1, 0.12),
                             greedy = FALSE) +
                          theme(plot.margin = margin(10, 0, 6, 15))

combined_4panel

ggsave(combined_4panel, 
       file = "Metabolites/Figures/Combined_Metabolites.png", width = 22, height = 12)