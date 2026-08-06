# Wearable Sensing Reveals the Structure of Cardiac Activation Associated with Everyday Driving

This repository contains the analysis code used to reproduce the results presented in the manuscript:

> **Wearable sensing reveals the structure of cardiac activation associated with everyday driving**  
> Yasmine Bouzid, MD Tanim Hasan, Michael Manser, and Ioannis Pavlidis

The repository accompanies the NUBI II naturalistic driving study and implements a baseline-referenced decomposition framework for understanding cardiac activation during everyday driving using continuous multimodal wearable sensing.

---

# Overview

Many physiologically meaningful behaviors are not characterized by rare extreme events but by modest physiological responses that occur repeatedly throughout daily life. Quantifying these repeated physiological exposures requires continuous measurements of physiology together with synchronized behavioral and environmental context.

Using one week of naturalistic monitoring from wearable devices, smartphones, and vehicles, this work characterizes the cardiac operating regime associated with everyday driving relative to non-driving sedentary behavior.

The proposed framework decomposes instantaneous heart rate into four interpretable components:

1. participant-specific physiological baseline,
2. context-specific cardiac offset,
3. additional modulation by behavioral, environmental, and individual characteristics,
4. residual variation.

The repository reproduces every analysis reported in the manuscript, including:

- exploratory data characterization,
- baseline-referenced operating-regime analysis,
- predictive decomposition using nested grouped cross-validation,
- elastic-net modulator analysis,
- long-horizon scenario projections,
- supplementary heart-rate missingness analyses,
- supplementary participant-level physiological traces.

---

# Repository structure

```text
.
├── Data/
│   └── NUBI_Data_60sec_Level_MASTER_CLEAN.csv
│
├── Scripts/
│   ├── 00_predictive_decomposition_nested_cv.R
│   ├── ...
│   ├── 07_figure7_table5_long_horizon_scenarios.R
│   ├── 08_supplementary_tableS1_hr_missingness_analysis.R
│   └── 09_supplementary_figureS1_participant_traces.R
│
├── Results/
│
├── README.md
├── LICENSE
└── DATA_USE.md
```

---

# Dataset

The repository includes the curated **60-second resolution** dataset used throughout the manuscript.

The dataset contains synchronized information from:

- Apple Watch heart rate
- Apple HealthKit physiological baseline
- smartphone sensing
- vehicle telemetry
- weather
- traffic
- psychometric instruments
- NASA-TLX workload measures

All analyses reported in the manuscript are reproduced from this dataset.

---

# Scientific workflow

The repository follows the analytical workflow of the manuscript.

```
Curated naturalistic dataset
            │
            ▼
Exploratory characterization
            │
            ▼
Baseline-referenced heart-rate transformation
            │
            ▼
Operating-regime analysis
            │
            ▼
Predictive decomposition
            │
            ▼
Elastic-net modulation analysis
            │
            ▼
Long-horizon scenario projections
            │
            ▼
Supplementary analyses
```

---

# Analysis pipeline

The principal scripts should be executed sequentially.

| Script | Primary output |
|---------|----------------|
| 00 | Nested grouped cross-validation and predictive decomposition |
| 01–06 | Figures 1–6 and associated manuscript tables |
| 07 | Figure 7 and long-horizon scenario analyses |
| 08 | Supplementary Table S1 (heart-rate missingness analyses) |
| 09 | Supplementary Figure S1 (participant-level HR traces) |

Each script creates its own output directory and reproduces the corresponding manuscript figures and tables.

---

# Methodological summary

The analysis uses participant-specific physiological baseline heart rate obtained from the Apple HealthKit `restingHeartRate` metric.

Baseline-referenced heart rate is defined as

```math
NHR = HR_{\mathrm{raw}} - HR_{\mathrm{base}}
```

The predictive framework represents instantaneous heart rate as

```math
HR_{\mathrm{raw}}
=
HR_{\mathrm{base}}
+
\tau_c
+
\phi(X)
+
\epsilon
```

where

- $HR_{\mathrm{base}}$ is the participant-specific baseline.
- $\tau_c$ is the context-specific offset.
- $\phi(X)$ represents modulation by observed covariates.
- $\epsilon$ is residual variation.

Prediction performance is evaluated using nested grouped cross-validation to ensure complete separation of participants between training and testing folds.

---

# Reproduced manuscript outputs

The repository reproduces all principal manuscript results, including

- Figure 1 – Exploratory sample characterization
- Figure 2 – Representative participant trace
- Figure 3 – Baseline-referenced operating regimes
- Figure 4 – Predictive decomposition framework
- Figure 5 – Predictive decomposition performance
- Figure 6 – Elastic-net modulators
- Figure 7 – Long-horizon scenario analyses

and

- Table 1 – Dataset inventory
- Table 2 – Missingness summary
- Table 3 – Operating-regime statistics
- Table 4 – Predictive decomposition performance
- Table 5 – Long-horizon scenario projections

The repository also reproduces

- Supplementary Table S1
- Supplementary Figure S1

---

## Software requirements

The analyses were developed in R. R version 4.6 or later is recommended.

The analysis scripts use the following R packages:

- `data.table`
- `doParallel`
- `dplyr`
- `forcats`
- `ggpattern`
- `ggplot2`
- `lubridate`
- `patchwork`
- `purrr`
- `readr`
- `rlang`
- `scales`
- `stringr`
- `tibble`
- `tidymodels`
- `tidyr`

Package dependencies are loaded by the individual scripts where required.

---

# Reproducibility

All scripts

- use deterministic random seeds where applicable,
- create output directories automatically,
- avoid manual intervention during execution,
- generate publication-ready figures directly from the released dataset.

Running the scripts sequentially reproduces the analyses reported in the accompanying manuscript.

---

# Data use

Please refer to **DATA_USE.md** for licensing terms, attribution requirements, and conditions governing reuse of the dataset.

---

# Citation

If you use this repository, please cite

Bouzid Y, Hasan MDT, Manser M, Pavlidis I.

*Wearable sensing reveals the structure of cardiac activation associated with everyday driving.*

---

# Acknowledgments

The NUBI II dataset was collected, curated, and prepared for public release by **MD Tanim Hasan** and **Ioannis Pavlidis**.

Repository organization, analysis scripts, and reproducible computational workflow were prepared by **Ioannis Pavlidis**.

---

# Contact

Ioannis Pavlidis

Affective and Data Computing Laboratory (ACDC)

University of Houston

Houston, Texas, USA

Email: ipavlidis@uh.edu
