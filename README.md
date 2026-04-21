### Exploring the effects of first year cold-dormancy on the phenotypic, molecular and metabolic responses in head-started Gopher tortoises (*Gopherus polyphemus*). ###

This repository contains data and code for exploring differences in: body size, growth rate, change in telomere lenght, and metabolic measures (mitochondrial DNA density, glucose, triglycerides and acetyl CoA) between two treatments of head-started Gopher tortoises from Alabama, US. The two treatments are: **Constant-Warmth**, i.e. animals that were raised under standard Alabama constant greenhouse conditions for the whole study period of 11 months, and **Cold-Dormancy**, i.e. animals that experienced a period of simulated winter dormancy for ~2 months. 
Data, code and results can be found within the folder for each respective analysis.

### Citation

Please use the following DOI to cite this data and associated analyses.

[![DOI](https://zenodo.org/badge/971060862.svg)](https://doi.org/10.5281/zenodo.15273407)




### Overview

```Gopher_tortoise_Head-start_2023/
├── Gopher_tortoise_Head-start_2023.Rproj
│
├── Growth_data/
│   ├── Code/
│   │   └── GT_Growth_Repeated_Measures.R    # Analysis for body size and growth rate differences
│   ├── Data/
│   │   ├── Growth_Metadata.xlsx             # Growth metadata (Tortoise IDs Nest IDs, Tanks, treatments, measurements, sampling timepoints)
│   │   └── Growth_data_Working.csv          # Cleaned growth data used in analysis
│   ├── Figures/                             # All the figures produced from the morphology analysis
│
├── Metabolites/
|   ├── Code/
│   │   └── GT_Metabolites_Analysis_Repeated_Measures.R    # Analysis for metabolite measures (glucose, triglycerides, glucose)
│   ├── Data/                                              # Masterdatasets for each metabolite separately + triglyceride dilution linearity
│   ├── Figures/                                           # All the figures produced from the metabolite measures
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
### Software

All the analyses were performed in Rstudio version 2024.12.1. The packages needed for each analysis are specified in the corresponding script.
