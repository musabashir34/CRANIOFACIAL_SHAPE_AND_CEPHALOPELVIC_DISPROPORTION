################################################################################
# confounder_controlled_shape_analysis.R
#
# Executes sequential confounder-controlled Procrustes ANOVA via procD.lm().
# Mathematically removes variation explained by maternal age, height, and parity
# to isolate the pure, unconfounded effect of CPD status on neonatal shape.
#
# REQUIRES: geomorph, dplyr, ggplot2
# install.packages(c("geomorph", "dplyr", "ggplot2"))
################################################################################

set.seed(42)
N_PERM <- 9999


## ======================== 1. ANTERIOR SHAPE CONTROL =========================

# Construct the specialized geomorph data frame for the anterior symmetric view
sex_numeric <- as.numeric(as.factor(as.character(master$SEX))) - 1 
ant_dataframe <- geomorph.data.frame(
  shape        = ant_sym_shape,       # Extracted object-symmetry shape array [15 x 2 x 160]
  group        = factor(ant_group),   # Clinical group factor ("S" or "C")
  maternal_age = as.numeric(master$`Maternal age`),
  maternal_ht  = as.numeric(master$`M Height(cm)`),
  parity       = as.numeric(master$Parity),
  sex          = sex_numeric  # Now passed safely as a plain 0/1 numeric vector
)

# CRITICAL LOGIC: Type I (Sequential) Sum of Squares formula.
# Listing maternal covariates FIRST forces the model to calculate and partition 
# away shape variation due to maternal traits BEFORE evaluating the 'group' factor.
ant_controlled_model <- procD.lm(
  shape ~ maternal_age + maternal_ht + parity + sex + group,
  data = ant_dataframe,
  iter = N_PERM,
  print.progress = FALSE
)


## ========================= 2. LATERAL SHAPE CONTROL =========================

# Construct the specialized geomorph data frame for the lateral profile view
lat_dataframe <- geomorph.data.frame(
  shape        = lat_shape,           # Aligned lateral coordinate shape array [14 x 2 x 160]
  group        = factor(lat_group),   # Clinical group factor ("S" or "C")
  maternal_age = as.numeric(master$`Maternal age`),
  maternal_ht  = as.numeric(master$`M Height(cm)`),
  parity       = as.numeric(master$Parity),
  sex          = sex_numeric  # Now passed safely as a plain 0/1 numeric vector
)

# Mirroring the strict sequential hierarchy for the lateral matrix
lat_controlled_model <- procD.lm(
  shape ~ maternal_age + maternal_ht + parity + sex + group,
  data = lat_dataframe,
  iter = N_PERM,
  print.progress = FALSE
)


## =================== 3. EXTRACTING ADJUSTED SHAPE RESIDUALS =================

# A. Flatten the raw shapes to 2D matrices to ensure absolute index safety
# This step automatically preserves unique specimen tags on the row names
ant_raw_2d <- two.d.array(ant_sym_shape)
lat_raw_2d <- two.d.array(lat_shape)

# B. Fit models restricted EXCLUSIVELY to the maternal confounders
ant_confounder_only <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex, data = ant_dataframe, iter = 0)
lat_confounder_only <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex, data = lat_dataframe, iter = 0)

# C. MATHEMATICAL EXTRACTION: Subtract the fitted maternal shapes from raw shapes
# This isolates the clean, unconfounded shape coordinates in 2D matrix space
ant_residuals_2d <- ant_raw_2d - ant_confounder_only$fitted
lat_residuals_2d <- lat_raw_2d - lat_confounder_only$fitted

# D. CRITICAL ARRAY STRUCTURING FOR PERFORMANCE METRICS SCRIPT:
# Fold the clean 2D matrices back into definitive 3D arrays (p x k x n).
# This preserves the full biological coordinate scale required by the grid engine.
ant_adjusted_residuals <- arrayspecs(
  ant_residuals_2d,
  p = dim(ant_sym_shape)[1], 
  k = dim(ant_sym_shape)[2]
)

lat_adjusted_residuals <- arrayspecs(
  lat_residuals_2d,
  p = dim(lat_shape)[1], 
  k = dim(lat_shape)[2]
)


## ============ 4. VISUALIZATION: UNCONFOUNDED SHAPE SCATTERPLOTS =============

# Run an orthogonal PCA directly on the maternal-clean 3D arrays
ant_clean_pca <- gm.prcomp(ant_adjusted_residuals)
lat_clean_pca <- gm.prcomp(lat_adjusted_residuals)

# Extract variance importance values for plot axes labelling
ant_pc_importance <- ant_clean_pca$importance$importance
lat_pc_importance <- lat_clean_pca$importance$importance

# Assemble plotting coordinates for ggplot2 mapping
plotting_df <- tibble(
  Ant_PC1 = ant_clean_pca$x[, 1],
  Ant_PC2 = ant_clean_pca$x[, 2],
  Lat_PC1 = lat_clean_pca$x[, 1],
  Lat_PC2 = lat_clean_pca$x[, 2],
  Group   = factor(ifelse(master$CPD == 1, "CPD (High Risk)", "Control (Vaginal)"))
)

# Plotting the Maternal-Clean Lateral Morphospace Separation
anterior_plot <- ggplot(plotting_df, aes(x = Ant_PC1, y = Ant_PC2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.68, linewidth = 1) + 
  scale_color_manual(values = c("darkorange", "deepskyblue4")) +
  labs(
    x = paste0("PC1 (", round(ant_pc_importance[1] * 100, 1), "% Variance)"),
    y = paste0("PC2 (", round(ant_pc_importance[2] * 100, 1), "% Variance)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

lateral_plot <- ggplot(plotting_df, aes(x = Lat_PC1, y = Lat_PC2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.68, linewidth = 1) + 
  scale_color_manual(values = c("darkorange", "deepskyblue4")) +
  labs(
    x = paste0("PC1 (", round(lat_pc_importance[1] * 100, 1), "% Variance)"),
    y = paste0("PC2 (", round(lat_pc_importance[2] * 100, 1), "% Variance)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

## ============ 5. PARAMETRIC MANCOVA USING 10 RETAINED PRINCIPAL COMPONENTS =============

# A. Extract the first 10 principal component scores

ant_retained_scores_for_mancova <- ant_pca$x[, 1:10]
lat_retained_scores_for_mancova <- lat_pca$x[, 1:10]

# B.  Reconstruct and Validate the Combined Parametric Data Tables

mancova_df_ant <- master %>%
  mutate(
    Maternal_Age    = as.numeric(`Maternal age`),
    Maternal_Height = as.numeric(`M Height(cm)`),
    Parity          = as.numeric(Parity),
    Sex             = factor(SEX, levels = c("Female", "Male")),
    Group           = factor(group, levels = c("C", "S")) # Clinical factor
  )

mancova_df_lat <- master %>%
  mutate(
    Maternal_Age    = as.numeric(`Maternal age`),
    Maternal_Height = as.numeric(`M Height(cm)`),
    Parity          = as.numeric(Parity),
    Sex             = factor(SEX, levels = c("Female", "Male")),
    Group           = factor(group, levels = c("C", "S")) # Clinical factor
  )

# C. Anterior Mancova
# Fit the sequential parametric MANCOVA model (Type I Sum of Squares)
ant_mancova_fit <- manova(
  ant_retained_scores_for_mancova ~ Maternal_Age + Maternal_Height + Parity + Sex + Group,
  data = mancova_df_ant
)
# Extract Pillai's Trace summary matrix
ant_mancova_summary <- summary(ant_mancova_fit, test = "Pillai")

# C. Lateral Mancova
# Fit the sequential parametric MANCOVA model (Type I Sum of Squares)
lat_mancova_fit <- manova(
  lat_retained_scores_for_mancova ~ Maternal_Age + Maternal_Height + Parity + Sex + Group,
  data = mancova_df_lat
)
# Extract Pillai's Trace summary matrix
lat_mancova_summary <- summary(lat_mancova_fit, test = "Pillai")

# Cache final statistics for dynamic R Markdown inline text integration
ant_group_pillai <- ant_mancova_summary$stats["Group", "Pillai"]
ant_group_f      <- ant_mancova_summary$stats["Group", "approx F"]
ant_group_p      <- ant_mancova_summary$stats["Group", "Pr(>F)"]

lat_group_pillai <- lat_mancova_summary$stats["Group", "Pillai"]
lat_group_f      <- lat_mancova_summary$stats["Group", "approx F"]
lat_group_p      <- lat_mancova_summary$stats["Group", "Pr(>F)"]

#====================================================================\n")
#   CANONICAL VARIATE ANALYSIS (CVA): DISTANCE EXTRACTION MODULE     \n")
#====================================================================\n\n")

## ===================== A. ANTERIOR VIEW CVA =====================


# Run the CVA on the retained anterior principal components
ant_cva_fit <- Morpho::CVA(
  ant_retained_scores_for_mancova, 
  groups = mancova_df_ant$Group, 
  rounds = 9999,      # High-resolution permutations to match your MANCOVA parameters
  cv = TRUE,          # FIXED: Explicitly enables cross-validation output calculation
  prior = NULL,
  plot = FALSE
)

# Extract distances and permutation p-values from the Morpho object
ant_mahal_dist  <- as.numeric(ant_cva_fit$Dist$GroupdistMaha)
ant_mahal_p     <- as.numeric(ant_cva_fit$Dist$probsMaha)
ant_procr_dist  <- as.numeric(ant_cva_fit$Dist$GroupdistEuclid)
ant_procr_p     <- as.numeric(ant_cva_fit$Dist$probsEuclid)


## ===================== B. LATERAL VIEW CVA =====================")

# Run the CVA on the retained lateral principal components
lat_cva_fit <- Morpho::CVA(
  lat_retained_scores_for_mancova, 
  groups = mancova_df_lat$Group, 
  rounds = 9999,
  cv = TRUE,          # FIXED: Explicitly enables cross-validation output calculation
  prior = NULL,
  plot = FALSE
)

# Extract lateral metrics
lat_mahal_dist  <- as.numeric(lat_cva_fit$Dist$GroupdistMaha)
lat_mahal_p     <- as.numeric(lat_cva_fit$Dist$probsMaha)
lat_procr_dist  <- as.numeric(lat_cva_fit$Dist$GroupdistEuclid)
lat_procr_p     <- as.numeric(lat_cva_fit$Dist$probsEuclid)

# Assemble the CVA distance table dynamically using the cached variables 
cva_table_data <- tibble(
  `Anatomical View` = c("**Anterior (Symmetric Facial Plane)**", "**Lateral (Neurocranial Profile View)**"),
  `Mahalanobis Distance` = c(round(ant_mahal_dist, 4), round(lat_mahal_dist, 4)),
  `Mahalanobis p-value` = c(
    ifelse(ant_mahal_p < 0.001, "< 0.001*", sprintf("%.4f*", ant_mahal_p)),
    ifelse(lat_mahal_p < 0.001, "< 0.001*", sprintf("%.4f*", lat_mahal_p))
  ),
  `Procrustes Distance` = c(round(ant_procr_dist, 4), round(lat_procr_dist, 4)),
  `Procrustes p-value` = c(
    ifelse(ant_procr_p < 0.001, "< 0.001*", sprintf("%.4f*", ant_procr_p)),
    ifelse(lat_procr_p < 0.001, "< 0.001*", sprintf("%.4f*", lat_procr_p))
  )
)

#  Extraction of the Leave-One-Out (Jackknife) Cross-Validated Classifications

# Pull the 160-element classification vectors natively from the CVA object slots
ant_cv_predictions <- ant_cva_fit$class
lat_cv_predictions <- lat_cva_fit$class

# 2. Pull the matching 160-element actual labels natively from the same object
ant_actual_groups  <- ant_cva_fit$groups
lat_actual_groups  <- lat_cva_fit$groups

# 3. Generate 2x2 Cross-Validated Confusion Matrices Natively
# Row 1/Col 1 is guaranteed to be "C", Row 2/Col 2 is guaranteed to be "S"
ant_cv_table <- table(Actual = ant_actual_groups, Predicted = ant_cv_predictions)
lat_cv_table <- table(Actual = lat_actual_groups, Predicted = lat_cv_predictions)

# 4. Standard Matrix Metric Extraction Engine via Fixed Coordinates
# [1,1] = True Control, [2,2] = True Case, [1,2] = False Case, [2,1] = False Control
calculate_cva_metrics_final = function(m_table) {
  m <- as.matrix(m_table)
  
  total_correct  <- sum(diag(m))
  total_samples  <- sum(m)
  cv_accuracy    <- (total_correct / total_samples) * 100
  
  # Row 2 is Actual 'S' (Cases), Col 2 is Predicted 'S'
  sens_pct <- (m[2, 2] / sum(m[2, ])) * 100
  
  # Row 1 is Actual 'C' (Controls), Col 1 is Predicted 'C'
  spec_pct <- (m[1, 1] / sum(m[1, ])) * 100
  
  list(total_correct = total_correct, total_samples = total_samples, 
       accuracy = cv_accuracy, sensitivity = sens_pct, specificity = spec_pct)
}

# Run execution summary passes safely
ant_metrics <- calculate_cva_metrics_final(ant_cv_table)
lat_metrics <- calculate_cva_metrics_final(lat_cv_table)

# 5. Format an Automated Tibble for Kable Presentation Printing
cva_accuracy_report_table <- tibble(
  `Anatomical View` = c("**Anterior (Symmetric Facial Plane)**", "**Lateral (Neurocranial Profile View)**"),
  `Total Samples (N)` = c(ant_metrics$total_samples, lat_metrics$total_samples),
  `Correct Classifications` = c(ant_metrics$total_correct, lat_metrics$total_correct),
  `Cross-Validated Accuracy (%)` = c(round(ant_metrics$accuracy, 1), round(lat_metrics$accuracy, 1)),
  `Sensitivity (True CPD, %)` = c(round(ant_metrics$sensitivity, 1), round(lat_metrics$sensitivity, 1)),
  `Specificity (True Control, %)` = c(round(ant_metrics$specificity, 1), round(lat_metrics$specificity, 1))
)

# Overlapping Cross-Validated Kernel Density distributions along the CV1 axis show the shape phenotypic segregation for both Anterior and Lateral profiles."}
# 1. Extract Cross-Validated CV1 Scores Natively from the Matrix
ant_cv_scores <- as.numeric(ant_cva_fit$CVcv[, 1])
lat_cv_scores <- as.numeric(lat_cva_fit$CVcv[, 1])

# 2. Build Tabular Data Frames for ggplot Ingestion
# Uses the exact actual groups factor pulled directly from the model object
plot_df_ant <- tibble(
  CV1   = ant_cv_scores,
  Group = factor(ant_actual_groups, levels = c("C", "S"), labels = c("Control (Vaginal)", "High Risk (CPD)"))
)

plot_df_lat <- tibble(
  CV1   = lat_cv_scores,
  Group = factor(lat_actual_groups, levels = c("C", "S"), labels = c("Control (Vaginal)", "High Risk (CPD)"))
)

# 3. Render Publication Density Plot: Anterior Face View
gg_cva_ant <- ggplot(plot_df_ant, aes(x = CV1, fill = Group, color = Group)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_rug(alpha = 0.7, sides = "b") + 
  scale_fill_manual(values = c("deepskyblue4", "darkorange")) +
  scale_color_manual(values = c("deepskyblue4", "darkorange")) +
  labs(
    title = "A. Anterior Facial Profile",
    x = "Canonical Axis 1 (CV1 Scores)",
    y = "Kernel Density Estimation"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

# 4. Render Publication Density Plot: Lateral Neurocranial Profile
gg_cva_lat <- ggplot(plot_df_lat, aes(x = CV1, fill = Group, color = Group)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_rug(alpha = 0.7, sides = "b") +
  scale_fill_manual(values = c("deepskyblue4", "darkorange")) +
  scale_color_manual(values = c("deepskyblue4", "darkorange")) +
  labs(
    title = "B. Lateral Neurocranial Profile ",
    x = "Canonical Axis 1 (CV1 Scores)",
    y = NULL 
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )







