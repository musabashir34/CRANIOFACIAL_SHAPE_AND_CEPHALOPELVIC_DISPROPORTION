# Reproducible Analysis Code

This repository contains the R Markdown source for:

> Bashir, M.A., Sabo, I. Craniofacial Shape and Size as Anatomical Correlates of Cephalopelvic Disproportion: A Geometric Morphometric Case-Control Study in Nigerian Neonates.

## Contents

- `Main_Manuscript.Rmd`: full analysis and manuscript source. Knits to a self-contained HTML report by default, and to a formatted Word (.docx) manuscript as a secondary format.
- `figures_source/`: static landmark-scheme photographs (Figures 1–2) referenced by the manuscript.
- Data-processing/analysis scripts sourced by the `.Rmd` (`data_ingestion.R`, `primary_gm_analysis_geomorph.R`, `build_master_dataset_geomorph.R`, `measurement_error.R`, `traditional_anthropometric_regression.R`, `confounder_controlled_shape_analysis.R`, `size_analysis.R`, `allometry_controlled_shape_analysis.R`, `transformation_grid.R`, `error_bubble_plot.R`) should sit alongside the `.Rmd` in this folder.
- `elsevier-harvard.csl`: the citation style language file used to format the references used in the article

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
rmarkdown::render("Main_Manuscript.Rmd")
```

This produces a single self-contained `Main_Manuscript.html` with all figures embedded (no external image files needed to view it).

**Word (.docx) manuscript:**

```r
rmarkdown::render("reproducible_cpd_article.Rmd", output_format = "officedown::rdocx_document")
```

Figure and table numbers, and all in-text cross-references, are generated automatically (via `bookdown`/`officedown` `\@ref()`) and stay consistent between both output formats.
