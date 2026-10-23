# Craniofacial Shape and Cephalopelvic Disproportion — Analysis Code and Data

Data, analysis code and R Markdown source for:

> Bashir MA, Sulaiman IS. Craniofacial Shape and Size as Anatomical Correlates of Cephalopelvic Disproportion: A Geometric Morphometric Case-Control Study in Nigerian Neonates.

The study compares neonatal craniofacial shape (two-dimensional geometric morphometrics on standardised anterior and lateral photographs) and conventional anthropometry between 80 neonates delivered by emergency caesarean section for cephalopelvic disproportion (CPD) and 80 uncomplicated vaginal deliveries at Jahun General Hospital, Jigawa State, Nigeria. All measurements were taken **after delivery**, so the results describe a postnatal association with CPD case status, not a predictive or screening tool.

## Repository contents

### Data

| File | Contents |
|---|---|
| `ANTERIOR__IMAGES_SAMPLE_AND_CONTROL.TPS` | 15 anterior landmarks × 160 neonates (TPSDig2) |
| `LATERAL_IMAGES_4_SAMPLE_AND_CONTROL.TPS` | 14 lateral landmarks × 160 neonates |
| `ANTERIOR_REPLICATE.TPS`, `LATERAL_REPLICATE.TPS` | Re-digitised subsample (16 per view) for measurement error |
| `Data_Sheet_of_Control.xlsx`, `Data_sheet_of_the_sample.xlsx` | Maternal and neonatal covariates for controls and CPD cases |
| `Data_Sheet_of_Control_Gestational_Age.csv`, `Data_Sheet_of_Sample_Gestational_Age.csv` | Gestational age (completed weeks), keyed by `S/N(ID)` |

**Specimen codes.** Each TPS `IMAGE=` name contains a code `[A|L][S|C]###`: view (A = anterior, L = lateral), group (S = CPD case, C = control) and a three-digit ID matching `S/N(ID)` in the covariate sheets.

**Units.** Landmark coordinates are in **pixels**; no scale bar was photographed, so centroid size is a relative, within-study measure. Birth weight is recorded in **pounds** (`Weight(ibs)`) and converted to kilograms in the code (1 lb = 0.45359237 kg). Head circumference is in cm; BPD and OFD in mm.

Photographs are not shared, to protect participants' privacy. Figures 1–2 (`landmark_scheme_anterior.jpeg`, `landmark_scheme_lateral.jpeg`) are shown with parental consent.

### Code

`Main_Manuscript.Rmd` sources the scripts below in this order. Each script relies on objects created by the ones before it.

| # | Script | What it does |
|---|---|---|
| 1 | `data_ingestion.R` | Reads both TPS files and covariate sheets; links them by specimen code; merges gestational age |
| 2 | `primary_gm_analysis_geomorph.R` | GPA (object-symmetry GPA for the anterior view, `bilat.symmetry`); centroid size; PCA |
| 3 | `build_master_dataset_geomorph.R` | Builds the linked master data frame |
| 4 | `measurement_error.R` | Two-way Procrustes ANOVA of original vs replicate digitisations; variance-component repeatability (ICC); per-landmark error |
| 5 | `traditional_anthropometric_regression.R` | Multivariable logistic regression (Table 1): Wald CIs, VIFs, Nagelkerke R² |
| 6 | `confounder_controlled_shape_analysis.R` | Covariate-adjusted Procrustes ANOVA; adjusted morphospaces (Figs S3–S4); MANOVA on 10 PCs; CVA distances and leave-one-out classification (Tables S4–S5, Fig S5) |
| 7 | `size_analysis.R` | Sequential ANCOVA of centroid size (Table S3) |
| 8 | `allometry_controlled_shape_analysis.R` | Procrustes ANOVA with log(centroid size) (Table 2); allometry plots (Fig 3) |
| 9 | `reviewer_response_analyses.R` | Group × log(CS) interaction; size-covariate sensitivity (Table S8); PC-retention sensitivity; allometric slope CIs; gestational age models (Tables S6–S7) |
| 10 | `landmark_exclusion_sensitivity.R` | Group effect with the most error-prone landmarks (zygion, opisthocranion) excluded |
| 11 | `transformation_grid.R` | Size- and covariate-adjusted TPS grids (Fig 4); exports the graphical abstract |
| 12 | `error_bubble_plot.R` | Per-landmark measurement-error maps (Figs S1–S2) |
| 13 | `participants_characteristics.R` | Participant characteristics by group (Table S1) |

All shape models use sequential (Type I) sums of squares with entry order maternal age, maternal height, parity, infant sex, [log centroid size], CPD group, so the group effect is conditional on the terms entered before it. Permutation tests use 9,999 iterations with `set.seed(42)`.

### Other files

- `mybibfile.bib`: references.
- `vancouver.csl`: numbered citation style used by *Morphologie*.

## Requirements

Analyses were run in **R 4.5.2** with **geomorph 4.1.0**, **RRPP 2.1.2** and **Morpho 2.13**.

```r
install.packages(c(
  "geomorph", "RRPP", "Morpho", "readxl", "dplyr", "tidyverse", "tidymodels",
  "glmnet", "pROC", "broom", "car", "vip", "gridExtra", "ggrepel",
  "flextable", "officer", "officedown", "bookdown", "rmarkdown", "knitr"
))
```

The YAML header builds `mybibfile.bib` from a Zotero library using [`rbbt`](https://github.com/paleolimbot/rbbt), which needs Zotero with Better BibTeX running. To knit without Zotero, change the header line to:

```yaml
bibliography: mybibfile.bib
```

## Reproducing the analysis

Keep all data files, scripts, `mybibfile.bib`, `vancouver.csl` and the two landmark images in the same folder as `Main_Manuscript.Rmd`, then:

```r
# HTML report (default)
rmarkdown::render("Main_Manuscript.Rmd")

# Word manuscript
rmarkdown::render("Main_Manuscript.Rmd", output_format = "officedown::rdocx_document")
```

The full run takes several minutes because of the permutation tests.

**Objects after knitting.** The Knit button runs in a separate R session, so analysis objects are not left in your console. To keep them for further work, render into the global environment instead:

```r
rmarkdown::render("Main_Manuscript.Rmd", envir = globalenv())
```

## Outputs

Knitting writes the following to the project folder:

| Output | Used for |
|---|---|
| `master_linked_dataset_geomorph.csv` | Linked analysis dataset |
| `table1_logistic_regression.csv` | Table 1 |
| `Supplementary_Table_01_Participant_Characteristics.csv` | Table S1 |
| `anterior_per_landmark_error_geomorph.csv`, `lateral_per_landmark_error_geomorph.csv` | Figs S1–S2 |
| `interaction_group_x_logsize.csv` | Slope-homogeneity test |
| `allometric_slope_cis.csv` | Per-group and common within-group allometric slopes |
| `sensitivity_alternative_size_covariates.csv` | Table S8 |
| `pc_variance_explained.csv`, `pc_retention_sensitivity.csv` | PC-retention sensitivity |
| `table_S6_logistic_with_gestational_age.csv` | Table S6 |
| `gestational_age_shape_models.csv` | Table S7 |
| `figures_publication/` | Submission figures (Figs 3–4, S1–S4) and graphical abstract |

## Licence

- **Code** (`.R`, `.Rmd`): [MIT License](LICENSE).
- **Data** (TPS, covariate and gestational-age files): [CC BY 4.0](LICENSE-DATA.md). Please cite the article when reusing them.
- **Landmark-scheme photographs**: not licensed for reuse; they show identifiable neonates and are published with parental consent for this study only.

## Contact

Musa Ado Bashir, Department of Anatomy, Bayero University Kano, Nigeria — musabashir34@gmail.com
