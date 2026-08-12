# Reproducible Analysis Code

This repository contains the R Markdown source for:

> Bashir, M. A., & Sabo, I. Postnatal Craniofacial Geometric Morphometrics and Traditional Anthropometry as Markers of Cephalopelvic Disproportion: A Comparative Case-Control Study. Submitted to The American Journal of Human Biology

## Contents

- `Main_Manuscript.Rmd` — full analysis and manuscript source. Knits to a self-contained HTML report by default, and to a formatted Word (.docx) manuscript as a secondary format.
- static landmark-scheme photographs (`landmark_scheme_anterior.jpeg` and `landmark_scheme_lateral.jpeg`) referenced by the manuscript.
- Digitized landmark data for anterior facial and lateral craniofacial views (`ANTERIOR__IMAGES_SAMPLE_AND_CONTROL.TPS` and `LATERAL_IMAGES_4_SAMPLE_AND_CONTROL.TPS`)
- repeat landmark data file for intra-observer digiting error assessment for both anterior facial and lateral craniofacial views (`ANTERIOR_REPLICATE.TPS` and `LATERAL_REPLICATE.TPS`)
- files recording neonatal anthropometric measurements and confounding maternal variables of age, parity and height for both control and CPD cases (`Data_Sheet_of_Control.xlsx` and `Data_sheet_of_the_sample.xlsx`)
- Data-processing/analysis scripts sourced by the `.Rmd` (`data_ingestion.R`, `primary_gm_analysis_geomorph.R`, `build_master_dataset_geomorph.R`, `measurement_error.R`, `traditional_anthropometric_regression.R`, `confounder_controlled_shape_analysis.R`, `size_analysis.R`, `allometry_controlled_shape_analysis.R`, `transformation_grid.R`, `error_bubble_plot.R`) should sit alongside the `.Rmd` in this folder.

## Requirements

R (≥ 4.2) with the following packages:

```r
install.packages(c(
  "geomorph", "readxl", "dplyr", "tidyverse", "tidymodels",
  "glmnet", "pROC", "broom", "car", "vip", "gridExtra", "ggrepel",
  "bookdown", "officedown", "officer", "flextable"
))
```

## Rendering

**HTML (default — for reading in a browser / GitHub Pages):**

```r
rmarkdown::render("reproducible_cpd_article.Rmd")
```

This produces a single self-contained `reproducible_cpd_article.html` with all figures embedded (no external image files needed to view it).

**Word (.docx) manuscript:**

```r
rmarkdown::render("reproducible_cpd_article.Rmd", output_format = "officedown::rdocx_document")
```

Figure and table numbers, and all in-text cross-references, are generated automatically (via `bookdown`/`officedown` `\@ref()`) and stay consistent between both output formats.
