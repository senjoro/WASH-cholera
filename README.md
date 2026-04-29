# WASH-cholera
Global Cholera Projections &amp; WASH Acceleration (2025-2054)

This repository contains the analytical framework and projection models for assessing the impact of Water, Sanitation, and Hygiene (WASH) interventions on global cholera risk under the SSP3-7.0 climate scenario.

## Overview
The study utilizes a **Negative Binomial Generalized Additive Model (GAM)** to:
1. Identify non-linear thresholds (e.g., the 80% WASH coverage effect).
2. Project cholera burden across 26 endemic countries under three scenarios: BAU, HTC, and INT.
3. Conduct **Target-seeking analysis** to identify the annual WASH acceleration required (e.g., 7%) to meet the GTFCC 2030 goals.

## Key Findings
- **Threshold Effect:** Protective gains from drinking water services accelerate significantly after crossing **80% coverage**.
- **Temporal Lag:** Intensive interventions (1% annual growth) show a temporary lag in efficacy before delivering ~50% reduction in the long term.
- **Acceleration Required:** A **7% annual absolute improvement** in WASH is necessary to achieve a 90% reduction in cases by 2030.

## Requirements
- R (>= 4.0)
- Packages: `mgcv`, `dplyr`, `ggplot2`, `data.table`, `scales`

## Usage
1. Load your historical dataset as `cholera_data`.
2. Run `vc_min.R` to train the GAM model and execute the 2030 target-seeking simulation.
3. The script will output the projected cases for various acceleration rates and visualize the pathway to elimination.

## Citation
If you use this code, please cite our manuscript:
> *Decoupling Global Cholera Risk from Climate Hazards: Non-linear Thresholds and Accelerated Infrastructure Pathways to 2054*
