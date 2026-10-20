# LCA Truleo Replication Repository

Replication code for latent class analysis of AZ DPS trooper attitudes toward body-worn camera auditing technology (Truleo).

## Overview

This repository contains R scripts to replicate the latent class analysis presented in a manuscript examining how trooper attitudes toward Truleo technology changed from pre-implementation to post-implementation. The analysis tests measurement invariance across waves and characterizes distinct officer attitude profiles.

## Directory Structure

```
lca-truleo-replication/
├── data/
│   └── pool.clean.rds              # Cleaned and merged survey data
├── scripts/
│   ├── 00_packages.R               # Load all required R packages
│   ├── 02_invariance_test.R        # Test measurement invariance (multi-group LCA)
│   ├── 03_wave_specific_lca.R      # Main analysis: fit independent models for Wave 1 and Wave 2
│   ├── 04_robustness_exclude_comprehension.R
│   ├── 05_robustness_collapse_likert.R
│   ├── 06_robustness_exclude_pilot.R
│   ├── 07_robustness_treatment_control.R
│   └── 08_robustness_treatment_control_troopers.R
└── output/                         # Generated figures and tables
    ├── main/                       # Primary manuscript results
    └── appendix/                   # Supplementary results
```

## Quick Start

1. **Set working directory** to `lca-truleo-replication/`

2. **Run scripts in order:**
   ```r
   source("scripts/00_packages.R")      # Load packages
   source("scripts/02_invariance_test.R")
   source("scripts/03_wave_specific_lca.R")
   source("scripts/04_robustness_exclude_comprehension.R")
   source("scripts/05_robustness_collapse_likert.R")
   source("scripts/06_robustness_exclude_pilot.R")
   source("scripts/07_robustness_treatment_control.R")
   source("scripts/08_robustness_treatment_control_troopers.R")
   ```

3. **Results** will be saved to `output/main/` and `output/appendix/` as figures (PNG) and tables (DOCX)

## Data

The analysis uses `data/pool.clean.rds`, a cleaned and merged dataset from Wave 1 (pre-implementation) and Wave 2 (post-implementation) surveys. This file is pre-processed and ready for analysis.

## Key Analyses

### Script 02: Measurement Invariance Test
Tests whether the latent class structure is invariant across waves using multi-group latent class analysis (glca package).

**Output:**
- Chi-square test comparing configural vs metric models
- RData file: `data/lca_invariance.RData`

### Script 03: Wave-Specific LCA (Main Analysis)
Fits independent four-class latent class models for each wave separately.

**Outputs:**
- Model enumeration tables (k=1 to k=5)
- Item-response probability figures
- Modal response tables
- Demographic stratification by class
- Item-level convergence analysis (t-tests, effect sizes)
- Convergence plot with overlay of t-test results

### Scripts 04-08: Robustness Checks
Verify that findings are robust to alternative specifications:
- **04**: Exclude comprehension items
- **05**: Collapse 5-point Likert to 3-point scale
- **06**: Exclude pilot study respondents
- **07**: Compare treatment vs control arms
- **08**: Compare treatment vs control (troopers only)

## Requirements

See `scripts/00_packages.R` for all package dependencies. Core packages:

- **tidyverse** — Data manipulation and visualization
- **poLCA** — Latent class analysis
- **glca** — Multi-group latent class analysis
- **effsize** — Effect size calculations (Cohen's d)
- **flextable** — Professional table formatting
- **rio** — Data import/export

## Output Files

All results are saved to `output/main/` and `output/appendix/`:


## Citation

If you use this code, please cite the original manuscript [citation information to be added].

## Questions

For questions about the analysis or code, please contact sbwatts@txstate.edu
