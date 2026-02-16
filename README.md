# Radiocarbon dating of *Persoonia arborea*

This repository contains R scripts and data to reproduce the analyses and figures from the manuscript "Dendroecology and radiocarbon show pulse recruitment after fire and persistence under logging in *Persoonia arborea*".

## R scripts

**tree_ring_vs_radiocarbon_age_analysis.R**
Compares reconstructed tree age from tree ring counts against radiocarbon dating ages to validate the reliability of tree ring reconstruction. Uses Bayesian measurement error regression to account for uncertainty in radiocarbon dating estimates.

**age_to_1.3m_height_analysis.R**
Estimates the time required for *Persoonia arborea* individuals to reach 1.3 m height following logging, based on regeneration data with known harvest dates.

**age_size_allometries.R**
Calibrates and compares multiple nonlinear growth models (power, Michaelis-Menten, Weibull, Gompertz, Chapman-Richards, Weibull, Bontemps-Duplat) to predict tree age from size measurements (DBH and height). Creates Fig X from the manuscript and Fig. Y from appendix.

**phenology_analysis.R**
Analyses the probability of flowering and fruiting as a function of tree size and season using generalised additive models. Creates Figs X, Y, Z from the manuscript and FIg. Y from appendix.

## Data files

**data_tree_ring_vs_radiocarbon.csv**
Individual trees measured using both tree ring counts and radiocarbon dating.
- `tree_ring_age_yr`: estimated age from tree ring count
- `f14_age_yr`: estimated age from radiocarbon dating
- `min_f14_age_yr`, `max_f14_age_yr`: radiocarbon dating age range
- `se_log_f14_age_yr`: standard error on log scale (for error-in-variable regression)

**data_age_to_1.3m_height.csv**
Height measurements of regenerating individuals following logging. Use to estimate the time required for individuals to reach 1.3 m height.
- `persoonia_regen_height_m`: measured height
- `max_tree_age_yrs`: maximum age based on harvest year

**data_age_dbh_allometry.csv**
Tree age versus diameter at breast height.
- `age_yr`: mean estimated age
- `se_age_yr`: standard error of age estimate
- `dbh_cm`: diameter at breast height
- `method`: age reconstruction method (regeneration, tree ring, or radiocarbon dating)

**data_age_height_allometry.csv**
Tree age versus height (same structure as DBH dataset but using `height_m` instead of `dbh_cm`).

**data_phenology.csv**
Flowering and fruiting observations from VicForests surveys. Each row is an individual observation.
- `flower_binary`: flowering status (0/1)
- `fruit_binary`: fruiting status (0/1)
- `year`: observation year
- `day_of_the_year`: Julian day
- `dbh_cm`: diameter at breast height
- `height_m`: tree height

## Software requirements

R version 4.1 or higher with packages:
- ggplot2
- mgcv
- brms
- cmdstanr
- minpack.lm
- gsl
- data.table

- viridis

