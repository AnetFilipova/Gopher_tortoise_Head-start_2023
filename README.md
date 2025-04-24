### Exploring the effects of first year cold-dormancy on the phenotypic, molecular and metabolic responses in head-started Gopher tortoises (*Gopherus polyphemus*). ###

This repository contains data and code for estimating differences in: growth rate, telomere lenght, mitochondrial DNA density, and metabolic measures (glucose, acetyl CoA, triglycerides) between two groups of head-started Gopher tortoises from Alabama, US. The two treatment groups are: **Constant-Warmth**, i.e. animals that were raised under standard Alabama constant greenhouse conditions for the whole study period of 11 months, and **Cold-Dormancy**, i.e. animals that experienced a period of simulated winter dormancy for ~2 months. 
Data, code and results for each separate analysis are contained in this project, but there are independent Rmarkdown documents for each analysis that can be found within the folder for the respective analysis.

### Citation

Please use the following DOI to cite this data and associated functions/analyses.

[![DOI](https://zenodo.org/badge/942427473.svg)](https://doi.org/10.5281/zenodo.14977202)

### Overview

```Gopher_tortoise_Head-start_2023/
├── Gopher_tortoise_Head-start_2023.Rproj
│
├── Growth_data/
│   ├── Code/
│   │   └── GT_Cohort3_Growth.R              # Analysis for growth rate and body size differences
│   ├── Data/
│   │   ├── Growth_Metadata.xlsx             # Growth metadata (Nest_ID, measurements, etc.)
│   │   └── Growth_data_Working.csv          # Cleaned growth data used in analysis
│   ├── Figures/                             # All the figures produced from growth analysis
│   ├── GT_Growth.Rmd                        # R Markdown for growth/body size analysis
│   ├── GT_Growth.html                       # Rendered HTML output
│   ├── GT_Growth.pdf                        # Rendered PDF output of the growth/body size analysis
│   └── GT_Growth.md                         # Markdown export of the Rmd file
│
├── Metabolomics/
│   └── [metabolomics analysis files/scripts]
│
├── Telo_mtDNA/
│   ├── qPCR_PlateReviewData_2_25_25.xlsx    # Plate review metadata for qPCR
│   └── Manuscript_Analysis/
│       ├── Code/
│       │   └── [GT_MtDNA_Telo.R]            # Analysis for telomere length and mtDNA density differences
│       ├── Data/
│       │   └── [telomere/mtDNA data files]  # qPCR exported files, trait metadata, etc.
│       ├── Figures/                         # All figures generated from telomere/mtDNA analyses
│       ├── Telo_mtDNA.Rmd                   # R Markdown for telomeres/mtDNA analysis
│       ├── Telo_mtDNA.html                  # Rendered HTML output
│       ├── Telo_mtDNA.pdf                   # Rendered PDF output of the telomeres/mtDNA analysis
│       └── Telo_mtDNA.md                    # Markdown export of the Rmd file

```
