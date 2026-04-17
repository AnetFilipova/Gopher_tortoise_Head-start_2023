library(dplyr)
library(readr)
library(stringr)

# Clear memory
rm(list = ls(all = TRUE))

# Directory containing all CSVs

data_dir <- "C:/Users/aliam/Desktop/GT_qPCR/Telo"

# List all plate CSV files
all_files <- list.files(data_dir, pattern = "Plate.*Quantification_Cq_Results\\.csv$", full.names = TRUE)
all_files

# Empty list to store results
all_results_list <- list()

# Loop through each file
for (file in all_files) {
  
  plate <- read_csv(file) %>% 
    rename_all(str_trim) %>%          # Trim column names
    filter(str_detect(Sample, "^STD")) # Keep standards only
  
  if (!"Cq Mean" %in% colnames(plate)) {
    warning(paste("Column 'Cq Mean' not found in file:", file))
    next  # Skip this file
  }
  
  # Compute summary for this plate
  cv_summary <- plate %>%
    group_by(Sample) %>%
    summarise(
      Mean_Cq = mean(Cq, na.rm = TRUE),
      SD_Cq = sd(Mean_Cq, na.rm = TRUE),,
      CV_percent = (SD_Cq / Mean_Cq) * 100,
      .groups = "drop"
    ) %>%
    mutate(
      Plate = basename(file),
      Type = case_when(
        str_detect(Plate, "Telomeres") ~ "Telo",
        str_detect(Plate, "Multiplex") ~ "Mito",
        TRUE ~ "Unknown"
      )
    )
  
  # Store results in the list
  all_results_list[[basename(file)]] <- cv_summary
}

# Combine all plate results into one data frame
all_results <- bind_rows(all_results_list)
all_results

# Overall summary by type and standard
overall_summary <- all_results %>%
  group_by(Type, Sample) %>%
  summarise(
    Overall_Mean_Cq = mean(Mean_Cq, na.rm = TRUE),
    Overall_SD_Cq = sd(Mean_Cq, na.rm = TRUE),
    Overall_CV_percent = (Overall_SD_Cq / Overall_Mean_Cq) * 100,
    .groups = "drop"
  )

# View results
print(overall_summary)
