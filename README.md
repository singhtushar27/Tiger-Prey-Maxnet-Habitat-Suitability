# Tiger Prey Resource Suitability Modelling

A reproducible R-based workflow for modelling combined tiger prey-resource environmental suitability using Maxnet, GIS and spatial environmental data.

## Overview

This project presents a reproducible R implementation of a MaxEnt-style habitat suitability modelling workflow to identify environmental conditions associated with recorded tiger prey-resource observations.

The workflow integrates climatic, vegetation, terrain and accessibility-related environmental variables and uses the `maxnet` package in R to generate a continuous prey-resource suitability surface.

The resulting suitability layer is intended to support subsequent landscape-level ecological analysis, habitat assessment and tiger corridor research.

### Methodological Note

The original habitat modelling analysis was conducted using standalone MaxEnt software. This repository presents a reproducible R-based implementation using the `maxnet` package to document and reproduce the modelling workflow.

## Objectives

- Prepare and harmonize environmental predictor layers.
- Reduce multicollinearity among continuous environmental variables.
- Incorporate categorical land-use/land-cover information.
- Prepare prey-resource occurrence and background data.
- Develop a Maxnet habitat suitability model.
- Evaluate model performance using random holdout and spatial cross-validation.
- Generate a continuous prey-resource suitability raster.

## Environmental Predictors

The final modelling dataset contained 14 environmental predictors:

- Bio2 – Mean Diurnal Range
- Bio3 – Isothermality
- Bio7 – Annual Temperature Range
- Bio9 – Mean Temperature of Driest Quarter
- Bio10 – Mean Temperature of Warmest Quarter
- Bio12 – Annual Precipitation
- Bio17 – Precipitation of Driest Quarter
- Bio18 – Precipitation of Warmest Quarter
- NDVI
- Slope
- Distance to river
- Distance to road
- Distance to railway
- LULC

Continuous variables were assessed using Pearson correlation, with highly correlated predictors removed using a correlation threshold of 0.70. LULC was retained separately as a categorical predictor.

## Occurrence Data

The occurrence dataset contains observations representing multiple prey and resource categories, including:

- Cattle
- Chital / Spotted Deer
- Goat
- Hanuman Langur
- Hare
- Nilgai
- Peafowl
- Porcupine
- Rhesus Macaque
- Sambar
- Wild Pig

After removing duplicate coordinates and retaining one occurrence per environmental raster cell, 185 unique presence cells were used for modelling.

A background sample of 10,000 environmentally complete cells was used for model calibration.

## Modelling

The model was developed using the R package `maxnet`.

The final model used:

- 185 presence cells
- 10,000 background cells
- 14 environmental predictors
- LULC treated as a categorical variable
- Maxnet feature classes: linear, quadratic, product and hinge
- Regularization multiplier: 1
- Cloglog output for continuous suitability prediction

The final prediction was generated as a continuous suitability surface ranging from 0 to 1.

## Model Evaluation

### Random Holdout Validation

- AUC: **0.9948**
- TSS: **0.9730**

### Five-Fold Spatial Cross-Validation

- Mean AUC: **0.9703 ± 0.0184**
- Mean TSS: **0.9044 ± 0.0577**

Spatial cross-validation was used to provide a more spatially rigorous assessment of model performance and reduce the influence of spatial dependence associated with a random data split.

## Results

The final model produced a continuous tiger prey-resource suitability surface.

The main visual outputs are provided in the `Figures/` directory, while numerical model evaluation results are available in the `results/` directory.

## Software and Methods

The workflow was developed using:

- R
- `terra`
- `maxnet`
- `pROC`
- GIS-based environmental raster processing
- MaxEnt-style species distribution modelling
- Spatial cross-validation

## Repository Structure

```text
Tiger-Prey-Maxnet-Habitat-Suitability/
│
├── Code/
│   └── Tiger_Prey_Maxnet_Reproducible_Workflow.R
│
├── Figures/
│   ├── Environmental_Variable_Correlation.png
│   └── Tiger_Prey_Resource_Habitat_Suitability.jpg
│
├── results/
│   ├── model_summary.csv
│   ├── spatial_validation_AUC.csv
│   └── spatial_validation_TSS.csv
│
└── README.md
