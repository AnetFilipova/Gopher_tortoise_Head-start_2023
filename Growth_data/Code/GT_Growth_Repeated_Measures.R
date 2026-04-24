##### Repeated Measures Analysis for Gopher Tortoise Growth Study #####

## Load necessary libraries ##
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)
library(car)

# Clear memory
rm(list=ls(all = TRUE))

# Define color palette
cbbPalette <- c("#0072B2", "#D55E00")

#### Data Preparation ####
datum <- read.csv("Growth_data/Data/Growth_data_Working.csv", na.strings = "na")

# Filter out non-viable individuals
datum <- subset(datum, !(Tortoise_ID %in% c("Not_Viable", "GT2023_N06.01","GT2023_N05.06", "GT2023_N15.04")))

##############################################################################
#### ANALYSIS: BODY SIZE (SCL) - Random Slopes Model ####
##############################################################################

# Prepare data in long format with experimental days
body_size_full <- datum %>%
  select(Tortoise_ID, Treatment, Nest_ID, Tank, 
         CL_10, CL_15, CL_16, CL_19, CL_20,
         ExpDay_measure10_Before, ExpDay_measure15_During, 
         ExpDay_measure16_3_Weeks_Post, ExpDay_measure19_3_Months_Post,
         ExpDay_measure20_1_Year_Post) %>%
  pivot_longer(
    cols = starts_with("CL_"), 
    names_to = "Timepoint", 
    values_to = "SCL"
  ) %>%
  mutate(
    Julian_day = case_when(
      Timepoint == "CL_10" ~ ExpDay_measure10_Before,
      Timepoint == "CL_15" ~ ExpDay_measure15_During,
      Timepoint == "CL_16" ~ ExpDay_measure16_3_Weeks_Post,
      Timepoint == "CL_19" ~ ExpDay_measure19_3_Months_Post,
      Timepoint == "CL_20" ~ ExpDay_measure20_1_Year_Post
    ),
    Julian_day_centered = Julian_day - min(Julian_day, na.rm = TRUE),
    Treatment = as.factor(Treatment),
    Nest_ID = as.factor(Nest_ID),
    Tank = as.factor(Tank),
    Tortoise_ID = as.factor(Tortoise_ID),
    Timepoint = factor(Timepoint, levels = c("CL_10", "CL_15", "CL_16", "CL_19", "CL_20"))
  ) %>%
  filter(!is.na(SCL) & !is.na(Julian_day)) %>%
  select(-starts_with("ExpDay_"))

#### Fit Random Slopes Model ####

# Scale Julian_day_centered to improve convergence
body_size_full <- body_size_full %>%
  mutate(Julian_day_scaled = scale(Julian_day_centered)[,1])

## Question: How fast are tortoises growing, i.e overall growth rates (slopes)?
## Do cold-dormancy and constant-warmth tortoises grow at different rates over time and if so by how much?

model_scl <- lmer(SCL ~ Treatment * Julian_day_scaled + 
             (1 | Nest_ID) +
             (Julian_day_scaled | Tortoise_ID), 
             data = body_size_full,
             control = lmerControl(optimizer = "bobyqa", 
             optCtrl = list(maxfun = 20000)))

# Model summary
summary(model_scl)
anova(model_scl, type = 3)

##############################################################################
#### Model Diagnostics - Check Assumptions ####
##############################################################################

# 1. NORMALITY OF RESIDUALS
par(mfrow = c(1, 2))
qqnorm(resid(model_scl), main = "Q-Q Plot: Residuals")
qqline(resid(model_scl), col = "red")
hist(resid(model_scl), main = "Histogram of Residuals", xlab = "Residuals", breaks = 20)
par(mfrow = c(1, 1))

# Shapiro-Wilk test for normality
shapiro.test(resid(model_scl))

plot(predict(model_scl),resid(model_scl))

# 2. HOMOGENEITY OF VARIANCE (Homoscedasticity)
plot(fitted(model_scl), resid(model_scl),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted (Check Homoscedasticity)")
abline(h = 0, col = "red", lwd = 2)

# Check variance by treatment
plot(body_size_full$Treatment, resid(model_scl),
     xlab = "Treatment", ylab = "Residuals",
     main = "Residuals by Treatment")

# Levene's test for homogeneity of variance
leveneTest(resid(model_scl) ~ Treatment, data = body_size_full)

# 3. OUTLIERS
# Standardized residuals
std_resid <- scale(resid(model_scl))
plot(std_resid, main = "Standardized Residuals", 
     ylab = "Standardized Residuals", xlab = "Observation")
abline(h = c(-3, 3), col = "red", lty = 2)

# Identify potential outliers (|std residual| > 3)
outliers <- which(abs(std_resid) > 3)
if(length(outliers) > 0) {
  print(paste("Potential outliers at observations:", paste(outliers, collapse = ", ")))
  print(body_size_full[outliers, c("Tortoise_ID", "Treatment", "Timepoint", "SCL")])
}

##############################################################################
#### Post-hoc Comparisons for Significant Interaction ####
##############################################################################

# 1. Compare growth rates (slopes) between treatments
emm_slopes <- emtrends(model_scl, ~ Treatment, var = "Julian_day_scaled")
print(emm_slopes)

pairs_slopes <- pairs(emm_slopes)
print(pairs_slopes)

# Convert scaled slopes back to original units (mm/day)
scale_factor <- sd(body_size_full$Julian_day_centered, na.rm = TRUE)
slopes_original <- summary(emm_slopes)
slopes_original$Julian_day_scaled.trend <- slopes_original$Julian_day_scaled.trend / scale_factor
slopes_original$SE <- slopes_original$SE / scale_factor
print(slopes_original)

## Question: At each specific measurement how different in size tortoises from the two treatments are?
## And does this difference change over time?
body_size_full <- body_size_full %>%
  mutate(Timepoint_numeric = case_when(
    Timepoint == "CL_10" ~ 0,
    Timepoint == "CL_15" ~ 1,
    Timepoint == "CL_16" ~ 2,
    Timepoint == "CL_19" ~ 3,
    Timepoint == "CL_20" ~ 4
  ))

## Timepoint_numeric is a continuous vs. Timepoint is categorical
model_scl_categorical <- lmer(SCL ~ Treatment * Timepoint + 
                         (1 | Nest_ID) + 
                         (Timepoint_numeric | Tortoise_ID), 
                         data = body_size_full,
                         control = lmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 20000)))

emm_by_time <- emmeans(model_scl_categorical, ~ Treatment | Timepoint)
print(emm_by_time)

pairs_by_time <- pairs(emm_by_time)
print(pairs_by_time)

# 3. Test effect of time within each treatment
emm_time_by_treatment <- emmeans(model_scl_categorical, ~ Timepoint | Treatment)
print(emm_time_by_treatment)

consec_comp <- contrast(emm_time_by_treatment, method = "consec")
print(consec_comp)

# 4. Calculate effect size
cd_rate <- slopes_original$Julian_day_scaled.trend[slopes_original$Treatment == "Cold dormancy"]
ch_rate <- slopes_original$Julian_day_scaled.trend[slopes_original$Treatment == "Constant heat"]
percent_diff <- ((ch_rate - cd_rate) / cd_rate) * 100

confint(pairs_by_time)
### The 95% C.I. for End-Dormancy is ± 4.95 ((-15.010-(-10.06)))
### The 95% C.I. for 3-Week Post-Dormancy is ± 6.34 (-11.62-(-17.963))
### The 95% C.I. for 3-Month Post-Dormancy is ± 7.93 (-7.41-(-15.339))
### The 95% C.I. for 1.5-Year Post-Dormancy is ± 11.17 ((-8.14-(-19.309)))

##############################################################################
#### Visualizations ####
##############################################################################


# Define color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Individual growth trajectories
growth_trajectories <- ggplot(body_size_full, aes(x = Julian_day_centered, y = SCL, 
                           group = Tortoise_ID, color = Treatment)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(group = Treatment), method = "lm", se = TRUE, linewidth = 1.5) +
  scale_color_manual(values = c("Cold dormancy" = cbbPalette[6],    # Blue #0072B2
                                "Constant heat" = cbbPalette[7]),    # Orange #D55E00
                     labels = c("Cold dormancy" = "Cold-Dormancy", 
                                "Constant heat" = "Constant-Warmth"),
                     name = "Treatment") +
  labs(x = "Days since start of experiment", 
       y = "Straight Carapace Length (mm)") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20, b = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.title.align = 0.5
  )

# Mean growth trajectories with significance
body_size_summary <- body_size_full %>%
  group_by(Treatment, Timepoint, Julian_day_centered) %>%
  summarise(
    mean_SCL = mean(SCL, na.rm = TRUE),
    se_SCL = sd(SCL, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
growth_trajectories

# Extract p-value for annotation
p_value <- summary(pairs_slopes)$p.value
sig_label <- ifelse(p_value < 0.0001, "****",
                    ifelse(p_value < 0.001, "***",
                           ifelse(p_value < 0.01, "**",
                                  ifelse(p_value < 0.05, "*", "ns"))))

# Boxplots by timepoint
body_size_plot <- body_size_full %>%
  mutate(
    Treatment_plot = ifelse(Timepoint == "CL_10", "Combined", as.character(Treatment)),
    Treatment_plot = factor(Treatment_plot, levels = c("Combined", "Cold dormancy", "Constant heat")),
    Timepoint_label = factor(Timepoint, 
                             levels = c("CL_10", "CL_15", "CL_16", "CL_19", "CL_20"),
                             labels = c("BD", "ED", "3wkPD", "3moPD", "1.5yrPD"))  # Updated labels
  )
# Updated colors to match growth rate plot (use cbbPalette[6] and [7])
colors_plot <- c("Combined" = "gray70", 
                 "Cold dormancy" = cbbPalette[6],  # Blue #0072B2
                 "Constant heat" = cbbPalette[7])   # Orange #D55E00

sig_by_timepoint <- summary(pairs_by_time)
sig_labels_box <- data.frame(
  Timepoint = c("CL_15", "CL_16", "CL_19", "CL_20"),
  Timepoint_label = factor(c("ED", "3wkPD", "3moPD", "1.5yrPD"),
                           levels = c("BD", "ED", "3wkPD", "3moPD", "1.5yrPD")),
  stars = sapply(2:5, function(i) {  # Changed from 1:4 to 2:5 to skip CL_10
    p <- sig_by_timepoint$p.value[i]
    ifelse(p < 0.0001, "****",
           ifelse(p < 0.001, "***",
                  ifelse(p < 0.01, "**",
                         ifelse(p < 0.05, "*", ""))))
  })
)

# Match the spacing from growth rate plot
max_y <- max(body_size_plot$SCL, na.rm = TRUE)
ypos <- max_y + 0.02 * max_y  # Asterisk position (2% above max)
line_y <- ypos - 0.005 * max_y  # Line position (0.5% below asterisks)

sig_labels_box$y_pos <- ypos
sig_labels_box$line_y <- line_y

# Add numeric x positions for lines (matching the timepoint positions)
sig_labels_box$x_num <- c(2, 3, 4, 5)  # Positions for Dormancy end, 3-Wk, 3-Mo, 1.5-Yr

body_size <- ggplot(body_size_plot, aes(x = Timepoint_label, y = SCL)) +
  # Before (CL_10) - gray filled box
  geom_boxplot(data = subset(body_size_plot, Timepoint == "CL_10"),
               aes(x = Timepoint_label, y = SCL),
               color = "black", fill = "gray70", outlier.shape = NA, 
               width = 0.6, linewidth = 0.6) +
  geom_point(data = subset(body_size_plot, Timepoint == "CL_10"),
             aes(x = Timepoint_label, y = SCL),
             position = position_jitter(width = 0.2), 
             size = 2, color = "black", fill = "white", shape = 21, stroke = 0.5) +
  
  # Other timepoints - white filled boxes with colored outlines
  geom_boxplot(data = subset(body_size_plot, Timepoint != "CL_10"),
               aes(color = Treatment_plot),
               fill = "white", outlier.shape = NA, width = 0.6, 
               position = position_dodge(0.8), linewidth = 0.6) +
  geom_jitter(data = subset(body_size_plot, Timepoint != "CL_10"),
              aes(color = Treatment_plot),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
              size = 2, alpha = 0.4) +
  
  # Significance lines (horizontal bars under asterisks)
  geom_segment(data = sig_labels_box,
               aes(x = x_num - 0.15, xend = x_num + 0.15, 
                   y = line_y, yend = line_y),
               color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  
  # Significance stars
  geom_text(data = sig_labels_box,
            aes(x = Timepoint_label, y = y_pos, label = stars),
            size = 7, fontface = "bold", inherit.aes = FALSE) +
  
  # Color scales (outlines only for non-Before timepoints)
  scale_color_manual(values = colors_plot,
                     labels = c("Cold dormancy" = "Cold-Dormancy", 
                                "Constant heat" = "Constant-Warmth"),
                     name = "Treatment") +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  
  labs(x = "Timepoint", 
       y = "Straight Carapace Length (mm)") +
  
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20, b = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.title.align = 0.5
  )

body_size

# Save file as PNG for final figure production
ggsave(growth_trajectories, file="Growth_data/Figures/Growth_Trajectories.png", width=9, height=7, dpi=600)
ggsave(body_size, file="Growth_data/Figures/Body_Size.png", width=9, height=7, dpi=600)

# Sample size by Nest_ID
body_size_full %>%
  distinct(Tortoise_ID, Nest_ID, Treatment) %>%
  count(Nest_ID, Treatment) %>%
print(n = Inf)

##############################################################################
#### ANALYSIS: GROWTH RATE across intervals ####
##############################################################################

# Prepare data in long format with growth rates
growth_rate_full <- datum %>%
  select(Tortoise_ID, Treatment, Nest_ID, Tank,
         Growth_rate_Before, Growth_rate_During, 
         Growth_rate_3_Weeks_Post, Growth_rate_3_Months_Post,
         Growth_rate_1.5_Year_Post) %>%
  pivot_longer(
    cols = starts_with("Growth_rate_"),
    names_to = "Interval",
    values_to = "Growth_rate"
  ) %>%
  mutate(
    Treatment = as.factor(Treatment),
    Nest_ID = as.factor(Nest_ID),
    Tank = as.factor(Tank),
    Tortoise_ID = as.factor(Tortoise_ID),
    Interval = factor(Interval, 
                      levels = c("Growth_rate_Before", "Growth_rate_During",
                                 "Growth_rate_3_Weeks_Post", "Growth_rate_3_Months_Post",
                                 "Growth_rate_1.5_Year_Post"))
  ) %>%
  filter(!is.na(Growth_rate))

## Question: Do treatments produce different growth rates within specific experimental intervals?
model_growth <- lmer(Growth_rate ~ Treatment * Interval + 
                 (1 | Tank) +
                 (1 | Nest_ID) + 
                 (1 | Tortoise_ID), 
                 data = growth_rate_full,
                 control = lmerControl(optimizer = "bobyqa"))

# Model summary
summary(model_growth)
anova(model_growth, type = 3)


##############################################################################
#### Model Diagnostics ####
##############################################################################

# 1. NORMALITY OF RESIDUALS
par(mfrow = c(1, 2))
qqnorm(resid(model_growth), main = "Q-Q Plot: Residuals")
qqline(resid(model_growth), col = "red")
hist(resid(model_growth), main = "Histogram of Residuals", xlab = "Residuals", breaks = 20)
par(mfrow = c(1, 1))

shapiro.test(resid(model_growth))

# 2. HOMOGENEITY OF VARIANCE
plot(fitted(model_growth), resid(model_growth),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lwd = 2)

plot(growth_rate_full$Treatment, resid(model_growth),
     xlab = "Treatment", ylab = "Residuals",
     main = "Residuals by Treatment")

leveneTest(resid(model_growth) ~ Treatment, data = growth_rate_full)

# 3. OUTLIERS
std_resid <- scale(resid(model_growth))
plot(std_resid, main = "Standardized Residuals", 
     ylab = "Standardized Residuals", xlab = "Observation")
abline(h = c(-3, 3), col = "red", lty = 2)

outliers <- which(abs(std_resid) > 3)
if(length(outliers) > 0) {
  print(paste("Potential outliers at observations:", paste(outliers, collapse = ", ")))
  print(growth_rate_full[outliers, c("Tortoise_ID", "Treatment", "Interval", "Growth_rate")])
}

##############################################################################
#### Post-hoc Comparisons ####
##############################################################################

# 1. Treatment differences within each interval
emm_by_interval <- emmeans(model_growth, ~ Treatment | Interval)
print(emm_by_interval)

pairs_by_interval <- pairs(emm_by_interval)
print(pairs_by_interval)

# 2. Interval differences within each treatment
emm_interval_by_treatment <- emmeans(model_growth, ~ Interval | Treatment)
print(emm_interval_by_treatment)

# Consecutive interval comparisons
consec_intervals <- contrast(emm_interval_by_treatment, method = "consec")
print(consec_intervals)

# 3. Overall treatment effect
emm_treatment <- emmeans(model_growth, ~ Treatment)
print(emm_treatment)
pairs(emm_treatment)

confint(pairs_by_interval)

### The 95% C.I. for Winter 2023 is ± 0.31 (-0.4062-(-0.717))
### The 95% C.I. for Early Sprin 2024 is ± 0.31 (-0.4865-(-0.797))
### The 95% C.I. for Late Spring 2024 is ± 0.31 (0.0832-(-0.228))
### ### The 95% C.I. for Summer 2025 is ± 0.48 (0.2068-(-0.271))

##############################################################################
#### Visualizations ####
##############################################################################

# Define color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##### Prepare the growth rate data in long format #####
growth_rate_long <- datum %>%
  filter(!Tortoise_ID %in% c("Not_Viable", "GT2023_N06.01", "GT2023_N05.06", "GT2023_N15.04")) %>%
  select(Tortoise_ID, Treatment, Nest_ID, Tank, 
         Growth_rate_Before, Growth_rate_During, Growth_rate_3_Weeks_Post, 
         Growth_rate_3_Months_Post, Growth_rate_1.5_Year_Post) %>%
  pivot_longer(
    cols = starts_with("Growth_rate_"),
    names_to = "Timepoint",
    values_to = "Growth_rate"
  ) %>%
  mutate(
    Timepoint = case_when(
      Timepoint == "Growth_rate_Before" ~ "Fall 2023\n(BD – Hatching)",
      Timepoint == "Growth_rate_During" ~ "Winter 2023\n(ED – BD)",
      Timepoint == "Growth_rate_3_Weeks_Post" ~ "Early Spring 2024\n(3wkPD – ED)",
      Timepoint == "Growth_rate_3_Months_Post" ~ "Late Spring 2024\n(3moPD – 3wkPD)",
      Timepoint == "Growth_rate_1.5_Year_Post" ~ "Summer 2025\n(1.5yrPD – 3moPD)",
      TRUE ~ Timepoint
    ),
    Timepoint = factor(
      Timepoint,
      levels = c(
        "Fall 2023\n(BD – Hatching)", 
        "Winter 2023\n(ED – BD)", 
        "Early Spring 2024\n(3wkPD – ED)", 
        "Late Spring 2024\n(3moPD – 3wkPD)", 
        "Summer 2025\n(1.5yrPD – 3moPD)"
      )
    )
  )

# Separate Fall 2023 from other timepoints
fall_2023 <- growth_rate_long %>% filter(Timepoint == "Fall 2023\n(BD – Hatching)")
other_timepoints <- growth_rate_long %>% filter(Timepoint != "Fall 2023\n(BD – Hatching)")

# Define y positions for significance
ypos <- max(growth_rate_long$Growth_rate, na.rm = TRUE) + 0.02
line_y <- ypos - 0.005

growth_rate <- ggplot() +
  # Fall 2023 - single gray box (combined)
  geom_boxplot(data = fall_2023, aes(x = Timepoint, y = Growth_rate),
               color = "black", fill = "gray70", outlier.shape = NA, 
               width = 0.6, linewidth = 0.6) +
  geom_point(data = fall_2023, aes(x = Timepoint, y = Growth_rate),
             position = position_jitter(width = 0.2), 
             size = 2, color = "black", fill = "white", shape = 21, stroke = 0.5) +
  
  # Other timepoints - separate boxes by Treatment (white fill, colored outlines)
  geom_boxplot(data = other_timepoints, aes(x = Timepoint, y = Growth_rate, color = Treatment),
               fill = "white", outlier.shape = NA, width = 0.6, 
               position = position_dodge(0.8), linewidth = 0.6) +
  geom_jitter(data = other_timepoints, aes(x = Timepoint, y = Growth_rate, color = Treatment),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
              size = 2, alpha = 0.4) +
  
  # Mean points and connecting lines for ALL timepoints
  stat_summary(data = growth_rate_long, 
               aes(x = Timepoint, y = Growth_rate, color = Treatment, group = Treatment),
               fun = mean, geom = "point", size = 2,
               position = position_dodge(width = 0.8)) +
  stat_summary(data = growth_rate_long, 
               aes(x = Timepoint, y = Growth_rate, color = Treatment, group = Treatment),
               fun = mean, geom = "line", linewidth = 0.6,
               position = position_dodge(width = 0.8)) +
  
  # Significance annotations
  # Winter 2023 (p < 0.0001)
  geom_segment(aes(x = 1.85, xend = 2.15, y = line_y, yend = line_y), 
               color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  annotate("text", x = 2, y = ypos, label = "****", size = 7, fontface = "bold") +
  # Early Spring 2024 (p = 0.0001)
  geom_segment(aes(x = 2.85, xend = 3.15, y = line_y, yend = line_y), 
               color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  annotate("text", x = 3, y = ypos, label = "****", size = 7, fontface = "bold") +
  
  # Theme and labels
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20, b = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.title.align = 0.5
  ) +
  
  # Colors: match body size plot style
  scale_color_manual(
    values = c("Cold dormancy" = cbbPalette[6], "Constant heat" = cbbPalette[7], "Combined" = "gray70"),
    name = "Treatment",
    breaks = c("Cold dormancy", "Constant heat", "Combined"),
    labels = c("Cold-Dormancy", "Constant-Warmth", "Combined")
  ) +
  
  labs(
    x = "Time Interval",
    y = "Growth rate (g/day)"
  )

### Show the plot ###
growth_rate

ggsave(growth_rate, file="Growth_data/Figures/Growth_Rate.png", width=9, height=7, dpi=600)

# Sample size by treatment
table(datum$Treatment)


###################################
# Test for treatment effect on survival (recapture in the wild) at 1.5 years post-dormancy
##################################

# Cold-Dormancy:   28 released, 9 recaptured, 19 not recaptured
# Constant-Warmth: 29 released, 11 recaptured, 18 not recaptured

survival_table <- matrix(c(9, 19,     # Cold-Dormancy row
                           11, 18),   # Constant-Warmth row
                         nrow = 2, 
                         byrow = TRUE,
                         dimnames = list(
                           Treatment = c("Cold-Dormancy", "Constant-Warmth"),
                           Status    = c("Recaptured", "Not recaptured")
                         ))

# View the table to confirm it looks right before testing
survival_table

## Fisher's exact test
# Tests whether recapture probability differs between treatments
# Returns: p-value, odds ratio, and 95% CI for the odds ratio
fisher.test(survival_table)

## Calculate the percentage of recaptured animals in the field

# Cold-Dormancy
9/28*100 ## 32.1%

# Constant-Warmth
11/29*100 ## %37.9%

######## Combine Body Size and Growth Rate into two-panel figure ########
library(patchwork)
library(cowplot)
library(grid)

# Remove legends from both panels
combined_body_growth <- (body_size + theme(legend.position = "none")) / 
  (growth_rate + theme(legend.position = "none")) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 25),
        plot.tag.position = c(0.02, 0.98),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 36, face = "bold"))

# Create a standalone legend by making a small dummy plot and extracting its legend
legend_plot <- ggplot(data.frame(
  Treatment = factor(c("Cold-Dormancy", "Constant-Warmth"), 
                     levels = c("Cold-Dormancy", "Constant-Warmth")),
  x = c(1, 2), y = c(1, 2)
), aes(x = Treatment, y = y, color = Treatment)) +
  geom_boxplot(fill = "white", linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Cold-Dormancy" = "#0072B2", "Constant-Warmth" = "#D55E00"),
    name = "Treatment"
  ) +
  theme(
    legend.title = element_text(size = 17, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key = element_blank(),
    legend.background = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank()
  )

# Extract just the legend from that dummy plot
shared_legend <- cowplot::get_legend(legend_plot)

# Combine: main figure on the left, legend on the right (wider legend column)
final_figure <- ggdraw() +
  draw_plot(combined_body_growth, x = 0, y = 0, width = 0.82, height = 1) +
  draw_plot(shared_legend, x = 0.82, y = 0.48, width = 0.18, height = 0.2)

final_figure

# Save
ggsave(final_figure, 
       file = "Growth_data/Figures/Combined_Body_Size_Growth_Rate.png",
       width = 14, height = 15, dpi = 600,
       bg = "white")
