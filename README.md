# Mugwort Research Repository

## Table of Contents

1. [Project Overview](#project-overview)
2. [Data](#data)
3. [Script](#script)
4. [Dependencies and Versions](#dependencies-and-versions)
5. [Reproducibility Guidelines](#reproducibility-guidelines)
6. [Contact](#contact)

## Project Overview

This repository contains the data and analysis code for a research project examining the interactions between mugwort (Artemisia), goldenrod (Solidago), and grasshoppers in experimental mesocosm and laboratory settings. The project investigates:

- **Elemental content**: Carbon and nitrogen content of plant tissues and soil
- **Plant structure**: Leaf area, stem complexity, and structural characteristics of Artemisia and Solidago
- **Mass-gain experiments**: Laboratory experiments examining grasshopper growth rates, molting rates, and survival on different vegetation diets (Artemisia, Solidago, and Grass)
- **Mesocosm survival**: Field mesocosm experiments tracking grasshopper survival across different vegetation communities and predator treatments
- **Biomass effects**: Analysis of plant biomass responses to herbivory and predator presence across different vegetation treatments
- **Habitat domains**: Behavioral observations of grasshopper height, habitat use, and foraging substrate preferences in response to vegetation and predator treatments

The analysis includes statistical modeling using linear models, mixed-effects models, survival analysis, and generalized linear mixed models, with comprehensive model diagnostics and post-hoc comparisons.

## Data

The `Data/` directory contains the following CSV files:

- **`HabitatDomains.csv`**: Behavioral observation data including grasshopper location (height, width), perch substrate, and behavior type across different treatment combinations
- **`InitialAllometry.csv`**: Initial allometric measurements of plants at the start of the experiment
- **`InitialCover.csv`**: Initial percent cover measurements for Artemisia, Solidago, and Grass by cage
- **`LeafScan.csv`**: Leaf area measurements from leaf scanning
- **`MassGain.csv`**: Laboratory mass-gain experiment data tracking individual grasshopper mass over time, including molting events and survival status
- **`mugwort_cn.csv`**: Carbon and nitrogen content data for plant tissues and soil samples
- **`PlantHarvest.csv`**: Final biomass measurements from plant harvest at the end of the experiment
- **`StemData.csv`**: Stem structural data including total leaf area, leaf number, stem intersections, and stem length
- **`Survival.csv`**: Mesocosm survival data tracking grasshopper counts over time with treatment information

All data files are in CSV format and can be read directly into R using `read.csv()`.

## Script

The main analysis script is **`mugwortscript.R`**, which contains all code necessary to reproduce the analyses and figures presented in the manuscript. The script is organized into the following sections:

1. **Elemental content analysis**: Linear models examining carbon and nitrogen content differences between species
2. **Plant structure analysis**: Analysis of stem complexity, leaf area, and stem intersections
3. **Mass-gain lab experiment**: Growth trajectory modeling, molting rate analysis, and survival analysis for laboratory experiments
4. **Mesocosm survival**: Beta GLMM analysis of interval survival rates in field mesocosms
5. **Biomass effects**: Linear models examining species-specific and total biomass responses to treatments
6. **Habitat domains**: Mixed-effects models and linear models examining grasshopper height, habitat use, and foraging substrate preferences

The script includes model diagnostics using DHARMa, post-hoc comparisons using emmeans, and figure generation using ggplot2. All figures are saved to the `Figures/` directory.

## Dependencies and Versions

This analysis requires R and the following R packages:

- **lme4** (>= 1.1-*): For fitting linear mixed-effects models
- **parameters** (>= *): For extracting model parameters
- **MuMIn** (>= *): For model selection and information-theoretic metrics
- **tidyverse** (>= 1.3.0): For data manipulation and visualization (includes dplyr, ggplot2, tidyr, readr, etc.)
- **ggthemes** (>= *): For additional ggplot2 themes
- **survival** (>= 3.2-*): For survival analysis
- **survminer** (>= *): For survival curve visualization
- **DHARMa** (>= 0.4.0): For model diagnostics using residual simulation
- **emmeans** (>= 1.7.0): For estimated marginal means and post-hoc comparisons
- **glmmTMB** (>= 1.0.0): For generalized linear mixed models with beta family
- **patchwork** (>= 1.1.0): For combining multiple ggplot2 figures


## Reproducibility

To reproduce the analyses and figures from this repository:

1. **Set up the working directory**: Ensure your R working directory is set to the repository root (the directory containing `mugwortscript.R` and the `Data/` folder).

2. **Install dependencies**: Install all required R packages as listed in the [Dependencies and Versions](#dependencies-and-versions) section.

3. **Load the script**: Open `mugwortscript.R` in R or RStudio.

4. **Run the analysis**: Execute the script from top to bottom. The script will:
   - Load all required packages
   - Read data from the `Data/` directory
   - Perform all statistical analyses
   - Generate diagnostic plots (displayed in R)
   - Create and save figures to the `Figures/` directory

5. **Figure output**: All figures are saved as PNG files in the `Figures/` directory with 300 DPI resolution. The script generates both individual figures and combined multi-panel figures.

**Important notes**:
- The script assumes data files are located in a `Data/` subdirectory relative to the working directory
- Figure outputs are saved to a `Figures/` subdirectory (this directory should exist or will be created automatically by `ggsave()`)
- Some model diagnostics are displayed interactively in R and may require manual inspection
- The script was last updated on 16 July 2025

## Contact

For questions about the code or data, please open an issue in this repository or contact Nathalie Sommer.
