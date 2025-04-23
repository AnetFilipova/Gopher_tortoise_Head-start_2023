##### Exploring the effect of first-year cold dormancy on growth trajectories in head-started Gopher tortoises #####


## Load necessary libraries ##
library(lme4) ## Will be used for mixed effect models with random and fixed effects 
library(emmeans) ## Will be used for post-hoc comparisons
library(ggplot2) ## Will be used for making graphs/data visualization 
library(dplyr) ## Will be used for specifying logical operators
library(ggpubr) ## Will be used for combining individual plots into one plot
library(tidyr) ## Will be used for reshaping the data for easier visualization later on
library(ggrepel) ## Will be used to label individual points on the plots

#Clear memory
rm(list=ls(all = TRUE)) 

#### Data Preparation ####

## Reading the .csv file and exploring the structure of the data ##
datum <- read.csv("Growth_data/Data/Growth_data_Working.csv", na.strings = "na")


head(datum)
str(datum)

## Filter out individuals where Tortoise_ID is equal to "Not_Viable" as well as specific IDs that need to be removed from the analyses ##
## Removing those individuals usin g the 'subset' function ##
## The logical operator %in% checks if elements of Tortoise_ID are present in the given vector c("Not_Viable", "GT2023_N05.03", "GT2023_N06.01", "GT2023_N05.06", "GT2023_N15.04") ##
## The operator ! returns the values from the column 'Tortoise_ID' that are NOT equal to any of the elements in the specified vector ##
datum <- subset(datum, !(Tortoise_ID %in% c("Not_Viable", "GT2023_N05.03", "GT2023_N06.01", "GT2023_N05.06", "GT2023_N15.04")))

# Confirm the filtered dataset
head(datum)

############### Data Analysis ##################

## Question 1: What is the growth rate of animals before cold dormancy treatment?
## Mixed effect model with treatment as a fixed effect and Nest_ID and Tank as a random effects ##
## If individuals from the same nest are more similar to each other than to individuals from different nests, treating Nest_ID as a random effect accounts for this correlation, avoiding pseudoreplication ##
## Statistical reason: If animals within the same tank grow similarly due to shared conditions, treating Tank as a random effect accounts for environmental clustering, preventing false inflation of significance ##

model_before <- lmer(Growth_rate_Before ~ Treatment + (1 | Nest_ID) + (1 | Tank), data = datum)

## Summarize the results
summary(model_before)

## Printing the 2.5% lower conficence limit and the 97.5% upper confidence limit to calculate the 95% confidence interval
confint(model_before) # The 95% CI is 0.05882025

## Post-hoc tests
emmeans(model_before, pairwise ~ Treatment)

#### Checking for outliers with z-scores > 3 meaning more than 3 standard deviations away from the mean, which is considered a stronger outlier ####
# Calculate the z-scores for the Growth_rate_Before variable
datum <- datum %>%
  mutate(z_score = (Growth_rate_Before - mean(Growth_rate_Before, na.rm = TRUE)) / sd(Growth_rate_Before, na.rm = TRUE))

# Filter out outliers based on Tortoise_ID
outliers_before_tortoise <- datum %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Growth_rate_Before, z_score)

# Print the outliers corresponding to individual Tortoise_ID
print(outliers_before_tortoise)

#### Plotting the data ####

## Definining an object for a colour-blind friendly pallette ##
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

before <- ggplot(datum, aes(x = Treatment, y = Growth_rate_Before, color = Treatment)) +  # Define aesthetics: x-axis as Treatment, y-axis as Growth_rate_Before, and color by Treatment
  geom_boxplot(position = position_dodge(0.85)) +  # Add boxplots with dodged positions to avoid overlap
  geom_jitter(width = 0.10, alpha = 0.5, size = 2) + # Add jittered points to show individual data points with some transparency defined by alpha = 0.5
  geom_text_repel(aes(label = Tank), size = 4, box.padding = 0.4, point.padding = 0.3, max.overlaps = 20) + # Add labels
  ylab("Growth Rate (g/day)") +  # Label the y-axis
  xlab("Treatment") +  # Label the x-axis
  scale_color_manual(values = c(cbbPalette[[6]], cbbPalette[[7]]), name = "", labels = c("", "")) + # Manually set colors and labels for the Treatment variable
  theme_classic() +  
  theme(strip.background = element_blank(), legend.position = "none", # Removing the default grey facet background as well as the legend by specifying 'none'
        axis.title = element_text(size = 18),  # Increase axis labels size
        axis.text = element_text(size = 14)    # Increase tick labels (treatment labels) size
  )  
before

# Save file as PNG for final figure production
ggsave(before, file="Growth_data/Growth_rate_Before.png", width=9, height=7, dpi=600)

####### Question 2: What is the effect of treatment on growth rate during dormancy? ########
# Mixed effect model with treatment as a fixed effect and Nest_ID and Tank as a random effects

model_during <- lmer(Growth_rate_During ~ Treatment + (1 | Nest_ID) + (1 | Tank) , data = datum)

# Summarize the results
summary(model_during)

## Printing the 2.5% lower conficence limit and the 97.5% upper confidence limit to calculate the 95% confidence interval
confint(model_during) # The 95% CI is 0.05889759

# Post-hoc tests
emmeans(model_during, pairwise ~ Treatment)


#### Checking for outliers with z-scores > 3 meaning more than 3 standard deviations away from the mean, which is considered a stronger outlier ####
# Calculate the z-scores for the Growth_rate_During variable
datum <- datum %>%
  mutate(z_score = (Growth_rate_During - mean(Growth_rate_During, na.rm = TRUE)) / sd(Growth_rate_During, na.rm = TRUE))

# Filter out outliers based on Tortoise_ID
outliers_during_tortoise <- datum %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Growth_rate_During, z_score)

# Print the outliers corresponding to individual Tortoise_ID
print(outliers_during_tortoise)


##### Plot the data #####
during <- ggplot(datum, aes(x = Treatment, y = Growth_rate_During, color = Treatment)) +  # Define aesthetics: x-axis as Treatment, y-axis as Growth_rate_During, and color by Treatment
  geom_boxplot(position = position_dodge(0.85)) +  # Add boxplots with dodged positions to avoid overlap
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 2) + # Add jittered points to show individual data points with some transparency defined by alpha = 0.5
  ylab("Growth Rate (g/day)") +  # Label the y-axis
  xlab("Treatment") +  # Label the x-axis
  scale_color_manual(values = c(cbbPalette[[6]], cbbPalette[[7]]), name = "", labels = c("", "")) + # Manually set colors and labels for the Treatment variable
  theme_classic() +  
  theme(strip.background = element_blank(), legend.position = "none", # Removing the default grey facet background as well as the legend by specifying 'none'
        axis.title = element_text(size = 18),  # Increase axis labels size
        axis.text = element_text(size = 14)    # Increase tick labels (treatment labels) size
  )  
during

# Save file as PNG for final figure production
ggsave(during, file="Growth_data/Growth_rate_During.png", width=9, height=7, dpi=600)

####### Question 3: What is the effect of treatment on growth rate 3 weeks post-cold dormancy? #######
model_3Weeks <- lmer(Growth_rate_3_Weeks_Post ~ Treatment + (1 | Nest_ID) + (1 | Tank) , data = datum)

# Summarize the results
summary(model_3Weeks)

# Post-hoc tests
emmeans(model_3Weeks, pairwise ~ Treatment)

#### Checking for outliers with z-scores > 3 meaning more than 3 standard deviations away from the mean, which is considered a stronger outlier ####
# Calculate the z-scores for the Growth_rate_3_Weeks_Post  variable
datum <- datum %>%
  mutate(z_score = (Growth_rate_3_Weeks_Post - mean(Growth_rate_3_Weeks_Post, na.rm = TRUE)) / sd(Growth_rate_3_Weeks_Post, na.rm = TRUE))

# Filter out outliers based on Tortoise_ID
outliers_3_Weeks_tortoise <- datum %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Growth_rate_3_Weeks_Post, z_score)

# Print the outliers corresponding to individual Tortoise_ID
print(outliers_3_Weeks_tortoise) ## Returns GT2023_N05.01 and GT2023_N17.03 as outliers

# Remove the outliers based on Tortoise_ID
datum_clean <- datum %>%
  filter(!Tortoise_ID %in% c("GT2023_N05.01", "GT2023_N17.03"))

# Re-run the mixed model without the outliers and without Nest_ID since it showed 0 variance
model_3Weeks_clean <- lmer(Growth_rate_3_Weeks_Post ~ Treatment  + (1 | Tank), data = datum_clean)

# Summarize the new model
summary(model_3Weeks_clean)

## Printing the 2.5% lower conficence limit and the 97.5% upper confidence limit to calculate the 95% confidence interval
confint(model_3Weeks_clean) ## The 95% CI is 0.2199788

# Post-hoc tests
emmeans(model_3Weeks_clean, pairwise ~ Treatment)


###### Plot the data ######

three_weeks <- ggplot(datum_clean, aes(x = Treatment, y = Growth_rate_3_Weeks_Post, color = Treatment)) +  # Define aesthetics: x-axis as Treatment, y-axis as Growth_rate_3_Weeks_Post, and color by Treatment
  geom_boxplot(position = position_dodge(0.85)) +  # Add boxplots with dodged positions to avoid overlap
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 2) + # Add jittered points to show individual data points with some transparency defined by alpha = 0.5
  ylab("Growth Rate (g/day)") +  # Label the y-axis
  xlab("Treatment") +  # Label the x-axis
  scale_color_manual(values = c(cbbPalette[[6]], cbbPalette[[7]]), name = "", labels = c("", "")) + # Manually set colors and labels for the Treatment variable
  theme_classic() +  
  theme(strip.background = element_blank(), legend.position = "none", # Removing the default grey facet background as well as the legend by specifying 'none'
        axis.title = element_text(size = 18),  # Increase axis labels size
        axis.text = element_text(size = 14)    # Increase tick labels (treatment labels) size
  )  
three_weeks

# Save file as PNG for final figure production
ggsave(three_weeks, file="Growth_data/Growth_rate_3Weeks.png", width=9, height=7, dpi=600)


######## Question 4: What is the effect of treatment on growth rate 3 months post-cold dormancy? ########
model_3Months <- lmer(Growth_rate_3_Months_Post ~ Treatment + (1 | Nest_ID) + (1| Tank), data = datum)

# Summarize the results
summary(model_3Months)

## Printing the 2.5% lower conficence limit and the 97.5% upper confidence limit to calculate the 95% confidence interval
confint(model_3Months) # The 95% CI is 0.1626287

# Post-hoc tests
emmeans(model_3Months, pairwise ~ Treatment)

#### Checking for outliers with z-scores > 3 meaning more than 3 standard deviations away from the mean, which is considered a stronger outlier ####
# Calculate the z-scores for the Growth_rate_3_Months_Post variable
datum <- datum %>%
  mutate(z_score = (Growth_rate_3_Months_Post - mean(Growth_rate_3_Months_Post, na.rm = TRUE)) / sd(Growth_rate_3_Months_Post, na.rm = TRUE))

# Filter out outliers based on Tortoise_ID
outliers_3_Months_tortoise <- datum %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Growth_rate_3_Months_Post, z_score)

# Print the outliers corresponding to individual Tortoise_ID
print(outliers_3_Months_tortoise)

###### Plot the data ########

three_months <- ggplot(datum, aes(x = Treatment, y = Growth_rate_3_Months_Post, color = Treatment)) +  # Define aesthetics: x-axis as Treatment, y-axis as Growth_rate_3_Weeks_Post, and color by Treatment
  geom_boxplot(position = position_dodge(0.85)) +  # Add boxplots with dodged positions to avoid overlap
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 2) + # Add jittered points to show individual data points with some transparency defined by alpha = 0.5
  ylab("Growth Rate (g/day)") +  # Label the y-axis
  xlab("Treatment") +  # Label the x-axis
  scale_color_manual(values = c(cbbPalette[[6]], cbbPalette[[7]]), name = "", labels = c("", "")) + # Manually set colors and labels for the Treatment variable
  theme_classic() +  
  theme(strip.background = element_blank(), legend.position = "none", # Removing the default grey facet background as well as the legend by specifying 'none'
        axis.title = element_text(size = 18),  # Increase axis labels size
        axis.text = element_text(size = 14)    # Increase tick labels (treatment labels) size
  )  
three_months

# Save file as PNG for final figure production
ggsave(three_months, file="Growth_data/Growth_rate_3Months.png", width=9, height=7, dpi=600)

######### Combine the results from all three models into one graph ##########

# Reshape data to long format
datum_long <- datum %>%
  pivot_longer(
    cols = c(Growth_rate_Before, Growth_rate_During, Growth_rate_3_Weeks_Post, Growth_rate_3_Months_Post),
    names_to = "Timepoint",
    values_to = "Growth_rate"
  )

# Adjust the Timepoint column for better labels
datum_long$Timepoint <- factor(
  datum_long$Timepoint,
  levels = c("Growth_rate_Before", "Growth_rate_During", "Growth_rate_3_Weeks_Post", "Growth_rate_3_Months_Post"),
  labels = c("Before", "During", "3-Wk.Post", "3-Mo.Post")
)


### Plot the combined data ###
combined <- ggplot(datum_long, aes(x = Timepoint, y = Growth_rate, color = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8), 
               fill = "white", alpha = 0.5, linewidth = 0.5) +  # Boxplot with dodging
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
              size = 2, alpha = 0.4) +  # Jitter points with dodging
  stat_summary(fun = mean, geom = "point", size = 2, aes(group = Treatment), 
               position = position_dodge(width = 0.8)) +  # Mean points, properly dodged
  stat_summary(fun = mean, geom = "line", aes(group = Treatment), linewidth = 0.6, 
               position = position_dodge(width = 0.8)) +  # Mean lines, properly dodged
  theme_classic () +
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
  scale_color_manual(
    values = cbbPalette[c(6, 7)], 
    name = "Treatment", 
    labels = c("Cold-Dormancy", "Constant-Warmth")
  ) +
  labs(
    x = "Timepoint",  
    y = "Growth rate (g/day)"
  )

combined

# Save file as PNG for final figure production
ggsave(combined, file = "Growth_data/Combined_Growth_rate.png", width = 12, height = 6, dpi = 600)



########## Body size across treatments and across timepoints #############

# Define list of tortoises to exclude
datum_filtered <- datum %>%
  filter(!Tortoise_ID %in% c("Not_Viable", "GT2023_N05.03", "GT2023_N06.01", "GT2023_N05.06", "GT2023_N15.04"))

# Reshape to long format and keep Nest_ID and Tank
datum_long <- datum_filtered %>%
  select(Tortoise_ID, Treatment, Nest_ID, Tank, CL_10, CL_15, CL_16, CL_19) %>%
  pivot_longer(cols = starts_with("CL_"), names_to = "Timepoint", values_to = "SCL") %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("CL_10", "CL_15", "CL_16", "CL_19")),
    Nest_ID = as.factor(Nest_ID),
    Tank = as.factor(Tank)
  )

####### Mixed-effect models for differences in body size (SCL) between treatments across the timepoints ###

##### Model Before #####
cl10_data <- datum_long %>% filter(Timepoint == "CL_10")

model_cl10 <- lmer(SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank), data = cl10_data)

summary(model_cl10)

emmeans_cl10 <- emmeans(model_cl10, pairwise ~ Treatment)

summary(emmeans_cl10$emmeans) ## Post-hoc comparison 
summary(emmeans_cl10$contrasts)

#### Checking for outliers in SCL ####

# Calculate the z-scores for the SCL (Straight Carapace Length) at CL_10 to test for outliers
datum_long_CL10 <- datum_long %>%
  filter(Timepoint == "CL_10")

datum_long_CL10 <- datum_long_CL10 %>%
  mutate(z_score = (SCL - mean(SCL, na.rm = TRUE)) / sd(SCL, na.rm = TRUE))

# Now filter the outliers for CL_10 specifically
outliers_SCL_CL10 <- datum_long_CL10 %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Timepoint, SCL, z_score)

print(outliers_SCL_CL10) ## there are no outliers

##### Model During #####
cl15_data <- datum_long %>% filter(Timepoint == "CL_15")

model_cl15 <- lmer(SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank), data = cl15_data)

summary(model_cl15)

emmeans_cl15 <- emmeans(model_cl15, pairwise ~ Treatment)

summary(emmeans_cl15$emmeans) ## Post-hoc comparison 
summary(emmeans_cl15$contrasts)

# Calculate the z-scores for the SCL at CL_15 to test for outliers
datum_long_CL15 <- datum_long %>%
  filter(Timepoint == "CL_15")

datum_long_CL15 <- datum_long_CL15 %>%
  mutate(z_score = (SCL - mean(SCL, na.rm = TRUE)) / sd(SCL, na.rm = TRUE))

# Now filter the outliers for CL_15 specifically
outliers_SCL_CL15 <- datum_long_CL15 %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Timepoint, SCL, z_score)

print(outliers_SCL_CL15) ## there are no outliers

##### Model 3 Weeks Post #####

cl16_data <- datum_long %>% filter(Timepoint == "CL_16")

model_cl16 <- lmer(SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank), data = cl16_data)

summary(model_cl16)

emmeans_cl16 <- emmeans(model_cl16, pairwise ~ Treatment)

summary(emmeans_cl16$emmeans) ## Post-hoc comparison 
summary(emmeans_cl16$contrasts)

# Calculate the z-scores for the SCL at CL_16 to test for outliers
datum_long_CL16 <- datum_long %>%
  filter(Timepoint == "CL_16")

datum_long_CL16 <- datum_long_CL16 %>%
  mutate(z_score = (SCL - mean(SCL, na.rm = TRUE)) / sd(SCL, na.rm = TRUE))

# Now filter the outliers for CL_16 specifically
outliers_SCL_CL16 <- datum_long_CL16 %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Timepoint, SCL, z_score)

print(outliers_SCL_CL16) ## there are no outliers

###### Model 3 Months Post #####
cl19_data <- datum_long %>% filter(Timepoint == "CL_19")

model_cl19 <- lmer(SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank), data = cl19_data)

summary(model_cl19)

emmeans_cl19 <- emmeans(model_cl19, pairwise ~ Treatment)

summary(emmeans_cl19$emmeans) ## Post-hoc comparison 
summary(emmeans_cl19$contrasts)

# Calculate the z-scores for the SCL at CL_19 to test for outliers
datum_long_CL19 <- datum_long %>%
  filter(Timepoint == "CL_19")

datum_long_CL19 <- datum_long_CL19 %>%
  mutate(z_score = (SCL - mean(SCL, na.rm = TRUE)) / sd(SCL, na.rm = TRUE))

# Now filter the outliers for CL_19 specifically
outliers_SCL_CL19 <- datum_long_CL19 %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Timepoint, SCL, z_score)

print(outliers_SCL_CL19) ## Tortoise GT2023_N05.01 was identified as an outlier so will be removed from the model


# Remove the outlier from the CL_19 dataset
datum_long_CL19_no_outlier <- datum_long_CL19 %>%
  filter(Tortoise_ID != "GT2023_N05.01")

# Rerun the model without the outlier
model_cl19_no_outlier <- lmer(SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank), data = datum_long_CL19_no_outlier)

# Summary of the model
summary(model_cl19_no_outlier)

# Post-hoc comparison (if needed)
emmeans_cl19_no_outlier <- emmeans(model_cl19_no_outlier, pairwise ~ Treatment)
summary(emmeans_cl19_no_outlier$emmeans)
summary(emmeans_cl19_no_outlier$contrasts)

### Plotting Body Size data #####

# Define color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Remove NAs and the specified individuals for plotting
datum_long_clean <- datum_long %>%
  filter(!is.na(SCL)) %>%  # Remove NAs in SCL
  filter(!(Tortoise_ID %in% c("Not_Viable", "GT2023_N05.03", "GT2023_N06.01", "GT2023_N05.06", "GT2023_N15.04", "GT2023_N05.01")))  # Remove specified individuals

# Summarize mean and standard error of SCL for plotting
summary_data <- datum_long_clean %>%
  group_by(Timepoint, Treatment) %>%
  summarise(
    mean_SCL = mean(SCL, na.rm = TRUE),
    se_SCL = sd(SCL, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Re-label timepoints for better axis display
timepoint_labels <- c("CL_10" = "Before", 
                      "CL_15" = "During", 
                      "CL_16" = "3-Wk.Post", 
                      "CL_19" = "3-Mo.Post")

# Create the matching body size boxplot
body_size_box <- ggplot(datum_long_clean, aes(x = Timepoint, y = SCL, color = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8),
               fill = "white", alpha = 0.5, linewidth = 0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
              size = 2, alpha = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 2, aes(group = Treatment),
               position = position_dodge(width = 0.8)) +
  stat_summary(fun = mean, geom = "line", aes(group = Treatment),
               linewidth = 0.6, position = position_dodge(width = 0.8)) +
  scale_color_manual(
    values = cbbPalette[c(6, 7)],
    name = "Treatment",
    labels = c("Cold-Dormancy", "Constant-Warmth")
  ) +
  scale_x_discrete(labels = timepoint_labels) +  # <- relabel the timepoints here
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
    y = "Mean SCL (mm)"
  )

body_size_box


# Save the plot as PNG for final figure production
ggsave(body_size_box, file = "Growth_data/Body_size.png", width = 12, height = 6, dpi = 600)

###### Combining figures for growth rate and body size together #######
final_figure <- ggarrange(
  body_size_box,  
  combined,   
  labels = c("A", "B"),  # Use capital letters for labels
  ncol = 2,  
  common.legend = TRUE,  
  legend = "right",  
  label.x = 0.2  # Move labels slightly to the right
)

final_figure

## Save as png Final Combine figure ##
ggsave(
  filename = "Growth_data/Final_Combined_Figure.png",
  plot = final_figure,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white",  # Fixes transparency issues
  device = "png"
)
