### Loading all the libraries that will be used in this code

``` r
library(lme4) ## Will be used for mixed effect models with random and fixed effects 
library(emmeans) ## Will be used for post-hoc comparisons
library(ggplot2) ## Will be used for making graphs/data visualization 
library(dplyr) ## Will be used for specifying logical operators
library(ggpubr) ## Will be used for combining individual plots into one plot
library(tidyr) ## Will be used for reshaping the data for easier visualization later on
library(ggrepel) ## Will be used to label individual points on the plots
```

# Data Preparation

``` r
datum <- read.csv("Growth_data_Working.csv", na.strings = "na") # Reading the data
```

### Filter out individuals where Tortoise_ID is equal to “Not_Viable” as well as specific IDs that need to be removed from the analyses. To remove those individuals, we will use the ‘subset’ function. The logical operator “%in%” checks if elements of Tortoise_ID are present in the given vector c(“Not_Viable”, “GT2023_N05.03”, “GT2023_N06.01”, “GT2023_N05.06”, “GT2023_N15.04”). The operator “!” returns the values from the column ‘Tortoise_ID’ that are NOT equal to any of the elements in the specified vector.

``` r
datum <- subset(datum, !(Tortoise_ID %in% c("Not_Viable", "GT2023_N05.03", 
"GT2023_N06.01", "GT2023_N05.06", "GT2023_N15.04")))
```

# Data Analysis

# Question 1: What is the growth rate of animals before cold dormancy treatment?

### We will run a linear mixed-effect model with Treatment as a fixed effect and Nest_ID and Tank as a random effects. If individuals from the same nest are more similar to each other than to individuals from different nests, treating Nest_ID as a random effect accounts for this correlation, avoiding pseudoreplication. The statistical reason is that if animals within the same tank grow similarly due to shared conditions, treating tank as a random effect accounts for environmental clustering, preventing false inflation of significance.

``` r
model_before <- lmer(Growth_rate_Before ~ Treatment + (1 | Nest_ID) +
    (1 | Tank), data = datum)

# Show summary statistics
summary(model_before)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: Growth_rate_Before ~ Treatment + (1 | Nest_ID) + (1 | Tank)
    ##    Data: datum
    ## 
    ## REML criterion at convergence: -132.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.78417 -0.65694  0.01712  0.65189  1.78394 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  Nest_ID  (Intercept) 0.0002578 0.01606 
    ##  Tank     (Intercept) 0.0016262 0.04033 
    ##  Residual             0.0034959 0.05913 
    ## Number of obs: 56, groups:  Nest_ID, 14; Tank, 10
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               0.38978    0.02173  17.937
    ## TreatmentConstant-Warmth -0.06098    0.03015  -2.023
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## TrtmntCns-W -0.693

``` r
# Printing the 2.5% lower confidence limit and the 97.5%
# upper confidence limit to calculate the 95% confidence
# interval.
confint(model_before)  # The 95% CI is 0.04529415
```

    ##                                2.5 %       97.5 %
    ## .sig01                    0.00000000  0.042705323
    ## .sig02                    0.01215073  0.066939834
    ## .sigma                    0.04772286  0.075287249
    ## (Intercept)               0.34745294  0.431946468
    ## TreatmentConstant-Warmth -0.11966254 -0.002022047

``` r
# Post-hoc comparisons between treatments
emmeans(model_before, pairwise ~ Treatment)
```

    ## $emmeans
    ##  Treatment       emmean     SE   df lower.CL upper.CL
    ##  Cold-Dormancy    0.390 0.0218 8.34    0.340    0.440
    ##  Constant-Warmth  0.329 0.0219 8.36    0.279    0.379
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ##  contrast                            estimate     SE   df t.ratio p.value
    ##  (Cold-Dormancy) - (Constant-Warmth)    0.061 0.0303 7.85   2.015  0.0794
    ## 
    ## Degrees-of-freedom method: kenward-roger

### Checking for outliers with z-scores \> 3 meaning more than 3 standard deviations away from the mean, which is considered a stronger outlier. This will be done for all the models below for each respective timepoint.

``` r
datum <- datum %>%
  mutate(z_score = (Growth_rate_Before - mean(Growth_rate_Before, na.rm = TRUE)) / 
  sd(Growth_rate_Before, na.rm = TRUE))

# Filter out outliers based on Tortoise_ID
outliers_before_tortoise <- datum %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Growth_rate_Before, z_score)

# Print the outliers corresponding to individual Tortoise_ID
print(outliers_before_tortoise) # There are no outliers
```

    ## [1] Tortoise_ID        Growth_rate_Before z_score           
    ## <0 rows> (or 0-length row.names)

# Plotting the data for Question 1

``` r
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
"#D55E00", "#CC79A7")

before <- ggplot(datum, aes(x = Treatment, y = Growth_rate_Before, color = Treatment)) +
  geom_boxplot(position = position_dodge(0.85)) + 
  geom_jitter(width = 0.10, alpha = 0.5, size = 2) +
  geom_text_repel(aes(label = Tank), size = 4, box.padding = 0.4, point.padding = 0.3, 
                  max.overlaps = 20) +
  ylab("Growth Rate (g/day)") +
  xlab("Treatment") +
  scale_color_manual(values = c(cbbPalette[[6]], cbbPalette[[7]]), name = "",
                     labels = c("", "")) +
  theme_classic() +  
  theme(strip.background = element_blank(), legend.position = "none",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)   
  )  
before
```

![](GT_Cohort3_Growth_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Save file as PNG for final figure production
ggsave(before, file="Growth_rate_Before.png", width=9, height=7, dpi=600)
```

# Question 2: What is the effect of treatment on growth rate during dormancy?

### Mixed effect model with Treatment as a fixed effect and Nest_ID and Tank as a random effects

``` r
model_during <- lmer(Growth_rate_During ~ Treatment + (1 | Nest_ID) + (1 | Tank),
data = datum)

summary(model_during)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: Growth_rate_During ~ Treatment + (1 | Nest_ID) + (1 | Tank)
    ##    Data: datum
    ## 
    ## REML criterion at convergence: -76.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5586 -0.3923 -0.0420  0.2397  3.6637 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  Nest_ID  (Intercept) 0.0031598 0.05621 
    ##  Tank     (Intercept) 0.0002478 0.01574 
    ##  Residual             0.0102632 0.10131 
    ## Number of obs: 56, groups:  Nest_ID, 14; Tank, 10
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               0.07822    0.02558   3.058
    ## TreatmentConstant-Warmth  0.54151    0.02931  18.474
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## TrtmntCns-W -0.568

``` r
confint(model_during) # The 95% CI is 0.05889759
```

    ##                               2.5 %    97.5 %
    ## .sig01                   0.00000000 0.1003494
    ## .sig02                   0.00000000 0.0584952
    ## .sigma                   0.08203378 0.1266309
    ## (Intercept)              0.02835180 0.1280811
    ## TreatmentConstant-Warmth 0.48167683 0.6000638

``` r
emmeans(model_during, pairwise ~ Treatment)
```

    ## $emmeans
    ##  Treatment       emmean     SE df lower.CL upper.CL
    ##  Cold-Dormancy   0.0782 0.0257 11   0.0216    0.135
    ##  Constant-Warmth 0.6197 0.0259 11   0.5627    0.677
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ##  contrast                            estimate     SE   df t.ratio p.value
    ##  (Cold-Dormancy) - (Constant-Warmth)   -0.542 0.0296 6.91 -18.287  <.0001
    ## 
    ## Degrees-of-freedom method: kenward-roger

``` r
datum <- datum %>%
  mutate(z_score = (Growth_rate_During - mean(Growth_rate_During, na.rm = TRUE)) /
  sd(Growth_rate_During, na.rm = TRUE))

outliers_during_tortoise <- datum %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Growth_rate_During, z_score)

print(outliers_during_tortoise) # There are no outliers
```

    ## [1] Tortoise_ID        Growth_rate_During z_score           
    ## <0 rows> (or 0-length row.names)

# Plotting the data for Question 2

``` r
during <- ggplot(datum, aes(x = Treatment, y = Growth_rate_During, color = Treatment)) +
  geom_boxplot(position = position_dodge(0.85)) +  
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 2) +
  ylab("Growth Rate (g/day)") +  # Label the y-axis
  xlab("Treatment") +  # Label the x-axis
  scale_color_manual(values = c(cbbPalette[[6]], cbbPalette[[7]]), name = "",
                     labels = c("", "")) +
  theme_classic() +  
  theme(strip.background = element_blank(), legend.position = "none",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)   
  )  
during
```

![](GT_Cohort3_Growth_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# Save file as PNG for final figure production
ggsave(during, file="Growth_rate_During.png", width=9, height=7, dpi=600)
```

# Question 3: What is the effect of treatment on growth rate 3 weeks post-cold dormancy?

### Mixed effect model with Treatment as a fixed effect and Nest_ID and Tank as a random effects

``` r
model_3Weeks <- lmer(Growth_rate_3_Weeks_Post ~ Treatment + (1 | Nest_ID) + (1 | Tank),
data = datum)

summary(model_3Weeks)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: Growth_rate_3_Weeks_Post ~ Treatment + (1 | Nest_ID) + (1 | Tank)
    ##    Data: datum
    ## 
    ## REML criterion at convergence: 46.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.5918 -0.4619 -0.1193  0.3216  3.3707 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Nest_ID  (Intercept) 0.01360  0.1166  
    ##  Tank     (Intercept) 0.03918  0.1979  
    ##  Residual             0.09193  0.3032  
    ## Number of obs: 56, groups:  Nest_ID, 14; Tank, 10
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)                0.6466     0.1105   5.851
    ## TreatmentConstant-Warmth   0.5915     0.1501   3.941
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## TrtmntCns-W -0.677

``` r
emmeans(model_3Weeks, pairwise ~ Treatment)
```

    ## $emmeans
    ##  Treatment       emmean    SE   df lower.CL upper.CL
    ##  Cold-Dormancy    0.647 0.111 8.85    0.395    0.898
    ##  Constant-Warmth  1.238 0.111 8.90    0.986    1.490
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ##  contrast                            estimate    SE   df t.ratio p.value
    ##  (Cold-Dormancy) - (Constant-Warmth)   -0.591 0.151 7.75  -3.925  0.0047
    ## 
    ## Degrees-of-freedom method: kenward-roger

``` r
datum <- datum %>%
  mutate(z_score = (Growth_rate_3_Weeks_Post - mean(Growth_rate_3_Weeks_Post, na.rm = TRUE)) /
  sd(Growth_rate_3_Weeks_Post, na.rm = TRUE))

outliers_3_Weeks_tortoise <- datum %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Growth_rate_3_Weeks_Post, z_score)

print(outliers_3_Weeks_tortoise)
```

    ##     Tortoise_ID Growth_rate_3_Weeks_Post  z_score
    ## 1 GT2023_N05.01                 2.521053 3.336697
    ## 2 GT2023_N17.03                 2.468421 3.225981

### There are two outliers, GT2023_N05.01 and GT2023_N17.03, so we need to filter them out and run the model again.

``` r
datum_clean <- datum %>%
  filter(!Tortoise_ID %in% c("GT2023_N05.01", "GT2023_N17.03"))

model_3Weeks_clean <- lmer(Growth_rate_3_Weeks_Post ~ Treatment  + (1 | Tank),
data = datum_clean)

summary(model_3Weeks_clean)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: Growth_rate_3_Weeks_Post ~ Treatment + (1 | Tank)
    ##    Data: datum_clean
    ## 
    ## REML criterion at convergence: 15.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.08279 -0.53538 -0.07705  0.37910  2.59762 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Tank     (Intercept) 0.02109  0.1452  
    ##  Residual             0.05946  0.2438  
    ## Number of obs: 54, groups:  Tank, 10
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               0.63853    0.07984   7.997
    ## TreatmentConstant-Warmth  0.51099    0.11363   4.497
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## TrtmntCns-W -0.703

``` r
confint(model_3Weeks_clean) ## The 95% CI is 0.2199788
```

    ##                              2.5 %    97.5 %
    ## .sig01                   0.0000000 0.2466127
    ## .sigma                   0.2005884 0.3051648
    ## (Intercept)              0.4840745 0.7934805
    ## TreatmentConstant-Warmth 0.2917524 0.7317100

``` r
emmeans(model_3Weeks_clean, pairwise ~ Treatment)
```

    ## $emmeans
    ##  Treatment       emmean     SE   df lower.CL upper.CL
    ##  Cold-Dormancy    0.639 0.0799 7.76    0.453    0.824
    ##  Constant-Warmth  1.150 0.0809 8.17    0.964    1.335
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ##  contrast                            estimate    SE   df t.ratio p.value
    ##  (Cold-Dormancy) - (Constant-Warmth)   -0.511 0.114 7.96  -4.493  0.0020
    ## 
    ## Degrees-of-freedom method: kenward-roger

# Plotting the data for Question 3

``` r
three_weeks <- ggplot(datum_clean, aes(x = Treatment, y = Growth_rate_3_Weeks_Post,
  color = Treatment)) +
  geom_boxplot(position = position_dodge(0.85)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 2) +
  ylab("Growth Rate (g/day)") +
  xlab("Treatment") +
  scale_color_manual(values = c(cbbPalette[[6]], cbbPalette[[7]]), name = "", labels = c("", "")) +
  theme_classic() +  
  theme(strip.background = element_blank(), legend.position = "none",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)
  )  
three_weeks
```

![](GT_Cohort3_Growth_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# Save file as PNG for final figure production
ggsave(three_weeks, file="Growth_rate_3Weeks.png", width=9, height=7, dpi=600)
```

# Question 4: What is the effect of treatment on growth rate 3 months post-cold dormancy?

### Mixed effect model with Treatment as a fixed effect and Nest_ID and Tank as a random effects

``` r
model_3Months <- lmer(Growth_rate_3_Months_Post ~ Treatment + (1 | Nest_ID) + (1| Tank),
data = datum)

summary(model_3Months)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: Growth_rate_3_Months_Post ~ Treatment + (1 | Nest_ID) + (1 |  
    ##     Tank)
    ##    Data: datum
    ## 
    ## REML criterion at convergence: 32.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.29961 -0.52573  0.07011  0.54660  2.14298 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Nest_ID  (Intercept) 0.016374 0.12796 
    ##  Tank     (Intercept) 0.001976 0.04446 
    ##  Residual             0.080598 0.28390 
    ## Number of obs: 56, groups:  Nest_ID, 14; Tank, 10
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               0.80451    0.06726   11.96
    ## TreatmentConstant-Warmth  0.01884    0.08200    0.23
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## TrtmntCns-W -0.606

``` r
confint(model_3Months) # The 95% CI is 0.1626287
```

    ##                               2.5 %    97.5 %
    ## .sig01                    0.0000000 0.2392906
    ## .sig02                    0.0000000 0.1566685
    ## .sigma                    0.2292185 0.3554776
    ## (Intercept)               0.6749433 0.9381987
    ## TreatmentConstant-Warmth -0.1441292 0.1811281

``` r
emmeans(model_3Months, pairwise ~ Treatment)
```

    ## $emmeans
    ##  Treatment       emmean     SE   df lower.CL upper.CL
    ##  Cold-Dormancy    0.805 0.0677 9.87    0.653    0.956
    ##  Constant-Warmth  0.823 0.0683 9.89    0.671    0.976
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ##  contrast                            estimate     SE   df t.ratio p.value
    ##  (Cold-Dormancy) - (Constant-Warmth)  -0.0188 0.0829 7.07  -0.227  0.8267
    ## 
    ## Degrees-of-freedom method: kenward-roger

``` r
datum <- datum %>%
  mutate(z_score = (Growth_rate_3_Months_Post - mean(Growth_rate_3_Months_Post, na.rm = TRUE)) /
  sd(Growth_rate_3_Months_Post, na.rm = TRUE))

outliers_3_Months_tortoise <- datum %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Growth_rate_3_Months_Post, z_score)

print(outliers_3_Months_tortoise)
```

    ## [1] Tortoise_ID               Growth_rate_3_Months_Post
    ## [3] z_score                  
    ## <0 rows> (or 0-length row.names)

# Plotting the data for Question 4

``` r
three_months <- ggplot(datum, aes(x = Treatment, y = Growth_rate_3_Months_Post,
  color = Treatment)) +
  geom_boxplot(position = position_dodge(0.85)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 2) +
  ylab("Growth Rate (g/day)") +
  xlab("Treatment") +
  scale_color_manual(values = c(cbbPalette[[6]], cbbPalette[[7]]), name = "",
  labels = c("", "")) +
  theme_classic() +  
  theme(strip.background = element_blank(), legend.position = "none",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)
  )  
three_months
```

![](GT_Cohort3_Growth_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
# Save file as PNG for final figure production
ggsave(three_months, file="Growth_rate_3Months.png", width=9, height=7, dpi=600)
```

# Now we want to combine the results from all four models into one graph. For this purpose, we will reshape the data to long format using the function pivot_longer.

``` r
datum_long <- datum %>%
  pivot_longer(
    cols = c(Growth_rate_Before, Growth_rate_During, Growth_rate_3_Weeks_Post,
    Growth_rate_3_Months_Post),
    names_to = "Timepoint",
    values_to = "Growth_rate"
  )

datum_long$Timepoint <- factor(
  datum_long$Timepoint,
  levels = c("Growth_rate_Before", "Growth_rate_During", "Growth_rate_3_Weeks_Post",
  "Growth_rate_3_Months_Post"),
  labels = c("Before", "During", "3-Wk.Post", "3-Mo.Post")
)
```

# Then we will plot the data from all four timepoints together and we will save it as a PNG for publication.

``` r
combined <- ggplot(datum_long, aes(x = Timepoint, y = Growth_rate, color = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8), 
               fill = "white", alpha = 0.5, linewidth = 0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
              size = 2, alpha = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 2, aes(group = Treatment), 
               position = position_dodge(width = 0.8)) +
  stat_summary(fun = mean, geom = "line", aes(group = Treatment), linewidth = 0.6, 
               position = position_dodge(width = 0.8)) +
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
```

![](GT_Cohort3_Growth_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
# Save file as PNG for final figure production
ggsave(combined, file = "Combined_Growth_rate.png", width = 12, height = 6, dpi = 600)
```

# The next thing we want to do is to estimate body size differences between the two treatments across the four timepoints. For this purpose, we will use the straight carapace length (mm) defined as “CL” in our data. The straight carapace length is a variable measured as straight line from the anterior (front) edge of the shell to the posterior (rear) edge, along the midline of the carapace. We will again filter out some individuals because they showed deviations from normal growth, hence they were excluded from the analyses.

``` r
datum_filtered <- datum %>%
  filter(!Tortoise_ID %in% c("Not_Viable", "GT2023_N05.03", "GT2023_N06.01",
  "GT2023_N05.06", "GT2023_N15.04"))

datum_long <- datum_filtered %>%
  select(Tortoise_ID, Treatment, Nest_ID, Tank, CL_10, CL_15, CL_16, CL_19) %>%
  pivot_longer(cols = starts_with("CL_"), names_to = "Timepoint", values_to = "SCL") %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("CL_10", "CL_15", "CL_16", "CL_19")),
    Nest_ID = as.factor(Nest_ID),
    Tank = as.factor(Tank)
  )
```

# Running the mixed-effect models for differences in body size (SCL) between treatments across the timepoints.

## Model for estimating body size when Carapace length is at measure 10. This corresponds to measure 10 which is the body size from hatching to ‘Before’, indicating before initiation of the cold dormancy treatment. At measurement 10 all animals were still under the same conditions.

``` r
cl10_data <- datum_long %>% filter(Timepoint == "CL_10")

model_cl10 <- lmer(SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank), data = cl10_data)

summary(model_cl10)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank)
    ##    Data: cl10_data
    ## 
    ## REML criterion at convergence: 257.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1352 -0.6067 -0.0028  0.5984  1.8068 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Nest_ID  (Intercept) 2.893    1.701   
    ##  Tank     (Intercept) 3.129    1.769   
    ##  Residual             3.309    1.819   
    ## Number of obs: 56, groups:  Nest_ID, 14; Tank, 10
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               65.6764     0.9792  67.074
    ## TreatmentConstant-Warmth  -2.9858     1.2279  -2.432
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## TrtmntCns-W -0.625

``` r
emmeans_cl10 <- emmeans(model_cl10, pairwise ~ Treatment)

summary(emmeans_cl10$emmeans)
```

    ##  Treatment       emmean    SE   df lower.CL upper.CL
    ##  Cold-Dormancy     65.7 0.980 11.5     63.5     67.8
    ##  Constant-Warmth   62.7 0.982 11.6     60.5     64.8
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95

``` r
summary(emmeans_cl10$contrasts)
```

    ##  contrast                            estimate   SE   df t.ratio p.value
    ##  (Cold-Dormancy) - (Constant-Warmth)     2.99 1.23 7.66   2.429  0.0426
    ## 
    ## Degrees-of-freedom method: kenward-roger

``` r
confint(model_cl10) # The 95% C.I. is 2.432148
```

    ##                               2.5 %     97.5 %
    ## .sig01                    0.9422280  2.8204918
    ## .sig02                    0.8702957  2.8682300
    ## .sigma                    1.4608881  2.3413323
    ## (Intercept)              63.7753542 67.5818112
    ## TreatmentConstant-Warmth -5.4246184 -0.5603227

``` r
datum_long_CL10 <- datum_long %>%
  filter(Timepoint == "CL_10")

datum_long_CL10 <- datum_long_CL10 %>%
  mutate(z_score = (SCL - mean(SCL, na.rm = TRUE)) / sd(SCL, na.rm = TRUE))

outliers_SCL_CL10 <- datum_long_CL10 %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Timepoint, SCL, z_score)

print(outliers_SCL_CL10) ## there are no outliers
```

    ## # A tibble: 0 × 4
    ## # ℹ 4 variables: Tortoise_ID <chr>, Timepoint <fct>, SCL <dbl>, z_score <dbl>

## Model for estimating body size when Carapace length is at measure 15. This corresponds to measure 15 which is the body size from ‘Before’ to ‘During’.’During’ indicates a measurement during the cold dormancy treatment but right at the end of the treatment. After that day the cold dormancy animals were prepared for arousal from dormancy.

``` r
cl15_data <- datum_long %>% filter(Timepoint == "CL_15")

model_cl15 <- lmer(SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank), data = cl15_data)

summary(model_cl15)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank)
    ##    Data: cl15_data
    ## 
    ## REML criterion at convergence: 306.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.9844 -0.6366 -0.1630  0.6433  2.7070 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Nest_ID  (Intercept)  4.126   2.0314  
    ##  Tank     (Intercept)  0.241   0.4909  
    ##  Residual             12.215   3.4950  
    ## Number of obs: 56, groups:  Nest_ID, 14; Tank, 10
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               69.7378     0.8915   78.23
    ## TreatmentConstant-Warmth  11.7298     1.0009   11.72
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## TrtmntCns-W -0.557

``` r
emmeans_cl15 <- emmeans(model_cl15, pairwise ~ Treatment)

summary(emmeans_cl15$emmeans)
```

    ##  Treatment       emmean    SE   df lower.CL upper.CL
    ##  Cold-Dormancy     69.7 0.896 11.3     67.8     71.7
    ##  Constant-Warmth   81.5 0.903 11.3     79.5     83.4
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95

``` r
summary(emmeans_cl15$contrasts)
```

    ##  contrast                            estimate   SE   df t.ratio p.value
    ##  (Cold-Dormancy) - (Constant-Warmth)    -11.7 1.01 6.86 -11.601  <.0001
    ## 
    ## Degrees-of-freedom method: kenward-roger

``` r
confint(model_cl15) ## The 95% C.I. is 2.016199
```

    ##                               2.5 %    97.5 %
    ## .sig01                    0.5973962  3.507058
    ## .sig02                    0.0000000  1.971003
    ## .sigma                    2.8288478  4.352135
    ## (Intercept)              68.0066897 71.487943
    ## TreatmentConstant-Warmth  9.6999580 13.732356

``` r
datum_long_CL15 <- datum_long %>%
  filter(Timepoint == "CL_15")

datum_long_CL15 <- datum_long_CL15 %>%
  mutate(z_score = (SCL - mean(SCL, na.rm = TRUE)) / sd(SCL, na.rm = TRUE))

outliers_SCL_CL15 <- datum_long_CL15 %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Timepoint, SCL, z_score)

print(outliers_SCL_CL15) ## there are no outliers
```

    ## # A tibble: 0 × 4
    ## # ℹ 4 variables: Tortoise_ID <chr>, Timepoint <fct>, SCL <dbl>, z_score <dbl>

## Model for estimating body size when Carapace length is at measure 16. This corresponds to measure 16 which is the body size from ‘During” to ’3 Weeks Post’. ‘3 Weeks Post’ indicates a measurement 3 weeks post- cold dormancy treatment. During this period the animals from the cold dormancy treatment are already living under the same conditions as the animals from the constant heat treatment

``` r
cl16_data <- datum_long %>% filter(Timepoint == "CL_16")

model_cl16 <- lmer(SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank), data = cl16_data)

summary(model_cl16)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank)
    ##    Data: cl16_data
    ## 
    ## REML criterion at convergence: 311
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.04910 -0.60698 -0.03726  0.62452  2.44199 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Nest_ID  (Intercept)  5.543   2.354   
    ##  Tank     (Intercept)  2.133   1.460   
    ##  Residual             11.618   3.408   
    ## Number of obs: 56, groups:  Nest_ID, 14; Tank, 10
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)                72.391      1.122   64.54
    ## TreatmentConstant-Warmth   13.880      1.315   10.55
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## TrtmntCns-W -0.583

``` r
emmeans_cl16 <- emmeans(model_cl16, pairwise ~ Treatment)

summary(emmeans_cl16$emmeans)
```

    ##  Treatment       emmean   SE   df lower.CL upper.CL
    ##  Cold-Dormancy     72.4 1.12 12.1     69.9     74.8
    ##  Constant-Warmth   86.3 1.13 12.2     83.8     88.7
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95

``` r
summary(emmeans_cl16$contrasts)
```

    ##  contrast                            estimate   SE   df t.ratio p.value
    ##  (Cold-Dormancy) - (Constant-Warmth)    -13.9 1.32 7.24 -10.500  <.0001
    ## 
    ## Degrees-of-freedom method: kenward-roger

``` r
confint(model_cl16) ## The 95% C.I. is 2.621092
```

    ##                              2.5 %    97.5 %
    ## .sig01                    1.031553  3.935275
    ## .sig02                    0.000000  2.850855
    ## .sigma                    2.750066  4.385224
    ## (Intercept)              70.220685 74.584305
    ## TreatmentConstant-Warmth 11.251507 16.493690

``` r
datum_long_CL16 <- datum_long %>%
  filter(Timepoint == "CL_16")

datum_long_CL16 <- datum_long_CL16 %>%
  mutate(z_score = (SCL - mean(SCL, na.rm = TRUE)) / sd(SCL, na.rm = TRUE))

outliers_SCL_CL16 <- datum_long_CL16 %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Timepoint, SCL, z_score)

print(outliers_SCL_CL16) ## there are no outliers
```

    ## # A tibble: 0 × 4
    ## # ℹ 4 variables: Tortoise_ID <chr>, Timepoint <fct>, SCL <dbl>, z_score <dbl>

## Model for estimating body size when Carapace length is at measure 19. This corresponds to measure 19 which is the body size from ‘3 Weeks Post’ to ‘3 Months Post’. ‘3 Months Post’ indicates a measurement 3 months post- cold dormancy treatment. This is also the last measurement taken for all animals before they were released in the wild.

``` r
cl19_data <- datum_long %>% filter(Timepoint == "CL_19")

model_cl19 <- lmer(SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank), data = cl19_data)

summary(model_cl19)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank)
    ##    Data: cl19_data
    ## 
    ## REML criterion at convergence: 360.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6402 -0.5761 -0.1328  0.6686  3.2450 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Nest_ID  (Intercept) 14.236   3.773   
    ##  Tank     (Intercept)  1.922   1.386   
    ##  Residual             30.870   5.556   
    ## Number of obs: 56, groups:  Nest_ID, 14; Tank, 10
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)                87.901      1.597  55.025
    ## TreatmentConstant-Warmth    9.398      1.754   5.357
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## TrtmntCns-W -0.545

``` r
emmeans_cl19 <- emmeans(model_cl19, pairwise ~ Treatment)

summary(emmeans_cl19$emmeans)
```

    ##  Treatment       emmean   SE   df lower.CL upper.CL
    ##  Cold-Dormancy     87.9 1.60 12.3     84.4     91.4
    ##  Constant-Warmth   97.3 1.61 12.4     93.8    100.8
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95

``` r
summary(emmeans_cl19$contrasts)
```

    ##  contrast                            estimate   SE   df t.ratio p.value
    ##  (Cold-Dormancy) - (Constant-Warmth)     -9.4 1.77 6.93  -5.315  0.0011
    ## 
    ## Degrees-of-freedom method: kenward-roger

``` r
datum_long_CL19 <- datum_long %>%
  filter(Timepoint == "CL_19")

datum_long_CL19 <- datum_long_CL19 %>%
  mutate(z_score = (SCL - mean(SCL, na.rm = TRUE)) / sd(SCL, na.rm = TRUE))

outliers_SCL_CL19 <- datum_long_CL19 %>%
  filter(abs(z_score) > 3) %>%
  select(Tortoise_ID, Timepoint, SCL, z_score)

print(outliers_SCL_CL19) ## Tortoise GT2023_N05.01 was identified as an outlier so will be removed from the model
```

    ## # A tibble: 1 × 4
    ##   Tortoise_ID   Timepoint   SCL z_score
    ##   <chr>         <fct>     <dbl>   <dbl>
    ## 1 GT2023_N05.01 CL_19      122.    3.47

``` r
# Remove the outlier from the CL_19 dataset
datum_long_CL19_no_outlier <- datum_long_CL19 %>%
  filter(Tortoise_ID != "GT2023_N05.01")

# Rerun the model without the outlier
model_cl19_no_outlier <- lmer(SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank),
data = datum_long_CL19_no_outlier)

summary(model_cl19_no_outlier)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: SCL ~ Treatment + (1 | Nest_ID) + (1 | Tank)
    ##    Data: datum_long_CL19_no_outlier
    ## 
    ## REML criterion at convergence: 339
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.86445 -0.61584 -0.09836  0.73417  2.29300 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Nest_ID  (Intercept) 10.078   3.175   
    ##  Tank     (Intercept)  4.564   2.136   
    ##  Residual             21.754   4.664   
    ## Number of obs: 55, groups:  Nest_ID, 14; Tank, 10
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)                87.929      1.565  56.192
    ## TreatmentConstant-Warmth    8.473      1.873   4.525
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## TrtmntCns-W -0.588

``` r
emmeans_cl19_no_outlier <- emmeans(model_cl19_no_outlier, pairwise ~ Treatment)
summary(emmeans_cl19_no_outlier$emmeans)
```

    ##  Treatment       emmean   SE   df lower.CL upper.CL
    ##  Cold-Dormancy     87.9 1.57 11.6     84.5     91.4
    ##  Constant-Warmth   96.4 1.59 12.1     92.9     99.9
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95

``` r
summary(emmeans_cl19_no_outlier$contrasts)
```

    ##  contrast                            estimate   SE  df t.ratio p.value
    ##  (Cold-Dormancy) - (Constant-Warmth)    -8.47 1.88 7.3  -4.502  0.0025
    ## 
    ## Degrees-of-freedom method: kenward-roger

``` r
confint(model_cl19_no_outlier) # The 95% C.I. is 3.730787
```

    ##                              2.5 %    97.5 %
    ## .sig01                    1.342439  5.366637
    ## .sig02                    0.000000  4.061572
    ## .sigma                    3.753833  6.004713
    ## (Intercept)              84.899469 90.997593
    ## TreatmentConstant-Warmth  4.722466 12.184039

# Plotting the data for body size across all four timepoints.

``` r
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
"#D55E00", "#CC79A7")

datum_long_clean <- datum_long %>%
  filter(!is.na(SCL)) %>%  # Remove NAs in SCL
  filter(!(Tortoise_ID %in% c("Not_Viable", "GT2023_N05.03", "GT2023_N06.01",
  "GT2023_N05.06", "GT2023_N15.04", "GT2023_N05.01")))  

summary_data <- datum_long_clean %>%
  group_by(Timepoint, Treatment) %>%
  summarise(
    mean_SCL = mean(SCL, na.rm = TRUE),
    se_SCL = sd(SCL, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

timepoint_labels <- c("CL_10" = "Before", 
                      "CL_15" = "During", 
                      "CL_16" = "3-Wk.Post", 
                      "CL_19" = "3-Mo.Post")

body_mass_box <- ggplot(datum_long_clean, aes(x = Timepoint, y = SCL, color = Treatment)) +
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
  scale_x_discrete(labels = timepoint_labels) +
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

body_mass_box
```

![](GT_Cohort3_Growth_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
# Save the plot as PNG for final figure production
ggsave(body_mass_box, file = "Body_mass.png", width = 12, height = 6, dpi = 600)
```

# Finally, we want to combine both figures for Growth rate and for Body size into one final figure with two panels. For this purpose, we will use the function ggarrange(). We will save the figure as a PNG for publication.

``` r
body_mass_box <- body_mass_box + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

combined <- combined + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

final_figure <- ggarrange(
  body_mass_box,  
  combined,   
  labels = c("A", "B"), 
  ncol = 2,  
  common.legend = TRUE,  
  legend = "right",  
  label.x = 0.05
)

final_figure
```

![](GT_Cohort3_Growth_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
ggsave(
  filename = "Final_Combined_Figure.png",
  plot = final_figure,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white", 
  device = "png"
)
```
