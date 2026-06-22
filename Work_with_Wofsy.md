# XCO2 / AOD Analysis — Script Reference

This directory contains R scripts for analyzing MethaneSAT (MSAT) XCO2 retrieval
biases using colocated VIIRS Deep Blue aerosol data. Scripts range from early
exploratory analyses to the final production pipeline.

---

## Core Pipelines

### `MSAT_VIIRS_model_pipeline.R`
The main unified pipeline. Colocates MethaneSAT and VIIRS at a 0.054° grid with
cosine-transformed viewing/solar angles, runs pooled AICc best-subset regression
(`dredge`) with scene ID as a fixed factor, applies VIF pruning for
multicollinearity, and outputs a full-bias-corrected XCO2 alongside per-scene
comparison maps and correction statistics.

### `MSAT_VIIRS_finalpipeline.R`
Advanced variant of the pipeline operating at finer 0.027° resolution. Adds
within-swath gap-filling via focal smoothing and Delaunay-triangulation
interpolation, across-track edge cropping (30 km buffer), and an interpolation
flag layer. Outputs a unified database, pooled dredge results, and per-scene
correction maps.

### `3.4.script.R`
Unified MSAT + VIIRS pipeline with cosine-transformed angles and full (non-scene-
mean-preserving) bias correction trained on pooled data. Uses a 0.05° grid,
within-swath interpolation via `akima`, and writes scaling parameter files for
consistent reapplication of the learned correction.

### `MSAT_VIIRS_whole_database.R`
Processes the entire MethaneSAT database with temporal VIIRS matching (±1.5 h).
Includes diagnostics on VIIRS/MSAT coverage percentages and overlap fractions,
uses 0.054° grid resolution, and logs all output to a timestamped file.

---

## Colocation & Correction Variants

### `MSAT_VIIRS_corrected.R`
Point-level colocation that preserves the native MSAT ~400 m pixel resolution by
extracting VIIRS values at each MSAT pixel location (±3 h temporal filter).
Applies a scene-mean-preserving XCO2 correction (zero global mean shift) and
outputs rasterized comparison plots at common resolution.

### `MSAT_VIIRS_corrected_costrans.R`
Variant of the above using cosine-transformed solar and viewing angles
(`cos(sza)`, `cos(vza)`, `cos(primary_angle)`) at the VIIRS 6 km native grid.
Produces publication-quality side-by-side maps of uncorrected and corrected XCO2.

### `msat_pipeline_debugged.R`
Debugging version of the processing pipeline with enhanced error reporting, better
NetCDF fill-value handling, relaxed minimum-pixel thresholds, and verbose
per-step logging. Produces the same corrected output while making failures easier
to diagnose.

---

## Statistical Modeling

### `XCO2_AOD_MLR_Workflow.R`
Integrated multi-level regression workflow using `dredge()` for AICc best-subset
selection with scene effects as fixed factors. Includes VIF-based predictor
pruning, a `speedglm` path for very large datasets, and PLS regression as a
comparative method. Outputs scene-level and pooled model coefficients.

### `AIC_XCO2_preds.R`
Analyzes XCO2 vs. delta pressure relationships stratified by aerosol type using
DP-weighted observations and prior shrinkage. Applies VIF-based predictor
selection and generates per-group regression coefficients with 95% confidence
intervals.

### `by_type_analysis.R`
Comprehensive stratified analysis by aerosol type. Methods include VIF pruning,
scene consistency F-tests, spatial autocorrelation testing with Moran's I
(KNN-based), spatial block cross-validation, and nonlinearity testing via quadratic
term comparisons with AIC model selection.

### `XCO2_correction.R`
General XCO2 correction framework applying aerosol-based bias removal using
colocated MSAT/VIIRS data. Rasterizes to VIIRS resolution, computes AOD at
multiple wavelengths, performs regression by aerosol type with multicollinearity
checks, and writes correction parameters.

### `XCO2_corr_coeff.R`
Correlation coefficient analysis examining relationships between XCO2 and
aerosol/atmospheric predictors by aerosol type. Computes AOD at custom wavelengths
via the Ångström power law, applies weighted regression, and reports per-group
correlation estimates with uncertainty.

### `XCO2_corr_norm.R`
Normalization-focused correction script that computes per-stratum delta-pressure
standard deviations for standardizing effect sizes across scenes and aerosol types.
Includes robust regression and VIF pruning.

---

## Bias & Aerosol Analysis

### `co2_ch4_bias_analysis_viirs_aod.R`
**Active script.** Analyzes retrieved CH4 and CO2 VCD biases against surface
albedo (1606 nm). For each MethaneSAT scene it: (1) locates or downloads the
colocated VIIRS Deep Blue granule within ±90 minutes, (2) extracts per-pixel AOD
(550 nm) and Ångström exponent, (3) computes per-aerosol-type mean AOD at 550 nm
and scaled to 1606 nm via the Ångström power law, and (4) annotates each output
plot (`bias_vs_albedo_<sid>.png`) with the top-3 aerosol types and their
AOD values in the bottom-right corner. GAM splines (`bs = "cs"`) are used for
the fitted bias curves.

### `co2_ch4_bias_analysis.R`
Baseline version of the above without VIIRS integration. Reads MethaneSAT L2
data, applies quality filtering, computes normalized percent bias in retrieved
CH4 and CO2 VCDs, and outputs GAM-smoothed bias-vs-albedo plots per scene.

### `bias_albedo_analysis.R`
Investigates XCO2 retrieval bias as a function of albedo at both the 1606 nm
and 1286 nm bands. Reads MSAT L2 data, computes normalized percent bias,
generates GAM-smoothed scatter plots, and annotates with top aerosol types
from colocated VIIRS data.

### `MSAT_VIIRS_analysis.R`
Early unified analysis script for MSAT/VIIRS colocation. Computes AOD at custom
wavelengths using the Ångström power law, defines aerosol type labels, and
performs regression by aerosol group with optional DP-weighted robust regression.

### `aerosol_comparison.R`
Compares aerosol optical properties across datasets or retrieval methods, examining
AOD and Ångström exponent consistency by aerosol type and evaluating systematic
differences in optical property estimates.

### `Global_plots.R`
Produces global-scale scatter plots of XCO2 vs. delta pressure and delta pressure
vs. AOD (1286 nm) with robust linear model (`rlm`) fitting. Uses non-parametric
bootstrapping (1000 iterations) for 95% confidence bands and density contours
to reveal clustering structure.

### `XCO2_vs_pred.R`
Comparative analysis of XCO2 against multiple predictors (delta pressure, albedo,
aerosol optical properties) stratified by aerosol type and Ångström exponent
category. Generates diagnostic plots with stratified regression overlays.

---

## Visualization

### `plotting_pipeline.R`
Visualization-focused script generating publication-quality output: side-by-side
XCO2 correction maps with consistent color scales, correction magnitude plots, and
cross-aerosol coefficient comparison figures with error bars and significance
indicators.

### `plotting_xco2_3.12.R`
Versioned plotting script (v3.12) producing global XCO2 maps, scatter plots with
fitted trends, and diagnostic figures with consistent color schemes and statistical
overlays (confidence intervals, density contours).

---

## Utilities

### `vars.R`
Shared configuration script. Defines directory paths, quality thresholds (AE
ranges, fill values), temporal matching windows, grid resolution parameters,
aerosol type labels, and other global constants used across the pipeline scripts.

### `test_quick.R`
Lightweight prototyping script for rapid iteration during development. Contains
simplified data extraction and minimal plotting for quickly validating pipeline
logic before full production runs.
