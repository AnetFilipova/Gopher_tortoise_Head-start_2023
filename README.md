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
│   ├── Figures/                                           # All the figures produced for the metabolite measures
│
├── Telo_mtDNA/
│   ├── Code/
│       ├── GT_MtDNA_Telo_Repeated_Measures.R              # Analysis for change in telomere length and mtDNA density from qPCR
│       └── GT_MtDNA_CV_Calculations.R                     # Analysis for calculating coefficient of variation (CV) between qPCR plates
│   ├── Data/                                              # qPCR exported files, trait metadata
│   ├── Figures/                                           # All figures generated for telomere/mtDNA

```
### Software

All the analyses were performed in Rstudio Version 2025.09.0. The packages needed for each analysis are specified in the corresponding script.
