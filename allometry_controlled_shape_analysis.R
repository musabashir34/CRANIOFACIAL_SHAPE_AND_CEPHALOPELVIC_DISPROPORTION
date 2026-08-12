################################################################################
# allometry_controlled_shape_analysis.R
#
# Tests if group shape differences persist after
# controlling for maternal confounders and infant sex AND Centroid Size (Allometry).
#
# REQUIRED PACKAGES: geomorph, dplyr, ggplot2
################################################################################

set.seed(42)
N_PERM <- 9999

#====================================================================\n")
#   ALLOMETRY & COVARIATE CONTROLLED PROCRUSTES ANOVA (procD.lm)     \n")
#====================================================================\n\n")


## ======================== 1. ANTERIOR VIEW ALLOMETRY ========================


# Reconstruct the geomorph data frame including log-transformed Centroid Size
ant_gdf_allom <- geomorph.data.frame(
  shape        = ant_sym_shape,
  group        = factor(ant_group),
  maternal_age = as.numeric(master$`Maternal age`),
  maternal_ht  = as.numeric(master$`M Height(cm)`),
  parity       = as.numeric(master$Parity),
  sex          = sex_numeric,
  log_size     = log(ant_csize) # Natural log of Centroid Size for standard allometric scaling
)

# STRICT HYERARCHICAL FORMULA: 
# Age + Height + Parity + Sex clears maternal  and infant gender variation. 
# log_size clears allometric scaling variation.
# group evaluates pure, size-independent shape differences.
ant_allom_model <- procD.lm(
  shape ~ maternal_age + maternal_ht + parity + sex + log_size + group,
  data = ant_gdf_allom,
  iter = N_PERM,
  print.progress = FALSE
)

## ======================== 2. LATERAL VIEW ALLOMETRY =========================


lat_gdf_allom <- geomorph.data.frame(
  shape        = lat_shape,
  group        = factor(lat_group),
  maternal_age = as.numeric(master$`Maternal age`),
  maternal_ht  = as.numeric(master$`M Height(cm)`),
  parity       = as.numeric(master$Parity),
  sex          = sex_numeric,
  log_size     = log(lat_csize)
)

lat_allom_model <- procD.lm(
  shape ~ maternal_age + maternal_ht + parity + sex + log_size + group,
  data = lat_gdf_allom,
  iter = N_PERM,
  print.progress = FALSE
)
ant_allom_summary <- as.data.frame(anova(ant_allom_model)$table)
lat_allom_summary <- as.data.frame(anova(lat_allom_model)$table)
# 2. Map tidy predictor rows to match your sequential design formulas perfectly
predictor_rows <- c("Maternal Age", "Maternal Height", "Parity", "Infant Sex", "log(Centroid Size)", "CPD Group Status", "Residuals")

# 3. Compile the Side-by-Side Publication Presentation Table
# Pulling precisely from the newly unpacked dataframe slots
side_by_side_anova_table <- tibble(
  `Predictor Source` = predictor_rows,
  
  # Anterior View Data Columns
  `Ant Df`    = ant_allom_summary$Df[1:7],
  `Ant SS`    = round(ant_allom_summary$SS[1:7], 5),
  `Ant F-value` = sapply(ant_allom_summary$F[1:7], format_f_value),
  `Ant p-value` = sapply(ant_allom_summary$`Pr(>F)`[1:7], format_p_value),
  
  # Lateral View Data Columns
  `Lat Df`    = lat_allom_summary$Df[1:7],
  `Lat SS`    = round(lat_allom_summary$SS[1:7], 5),
  `Lat F-value` = sapply(lat_allom_summary$F[1:7], format_f_value),
  `Lat p-value` = sapply(lat_allom_summary$`Pr(>F)`[1:7], format_p_value)
)



## ================= B. EXTRACTING PURE SIZE-CLEAN RESIDUALS ==================

# Fit models restricted EXCLUSIVELY to the maternal background confounders and log_size
ant_covar_only <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + log_size, data = ant_gdf_allom, iter = 0)
lat_covar_only <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + log_size, data = lat_gdf_allom, iter = 0)

# DEFINITIVE RESOLUTION: Explicitly cast the model's fitted lists into 2D matrices.
# This guarantees that the matrix layouts (Specimens x Coordinates) match perfectly,
# resolving the dimension errors during subtraction.
ant_fitted_matrix <- as.matrix(ant_confounder_only$fitted)
lat_fitted_matrix <- as.matrix(lat_confounder_only$fitted)

# Perform matrix subtraction to isolate the pure size-and-confounder clean shapes
ant_pure_clean_2d <- two.d.array(ant_sym_shape) - ant_fitted_matrix
lat_pure_clean_2d <- two.d.array(lat_shape) - lat_fitted_matrix

# Fold the clean 2D coordinates back into a standard 3D array (p x k x n)
# This format is required by the downstream manual RegScore and machine learning engines
ant_pure_residuals <- arrayspecs(ant_pure_clean_2d, p = dim(ant_sym_shape)[1], k = dim(ant_sym_shape)[2])
lat_pure_residuals <- arrayspecs(lat_pure_clean_2d, p = dim(lat_shape)[1], k = dim(lat_shape)[2])



## ====================== 5. VISUALIZING ALLOMETRIC SLOPES ====================

# --- RE-ENGINEERED MANUAL REGRESSION SCORE LOGIC (Drake & Klingenberg, 2008) ---
# To eliminate 'not a graphical parameter' warnings completely, we bypass plot()
# and compute scores natively via direct matrix multiplication.
# RegScore = Flattened_Coordinates %*% (Beta_Slope_Vector / ||Beta_Slope_Vector||)

compute_silent_reg_score <- function(model_obj, raw_shape_array) {
  # 1. Convert the 3D landmark array into a flattened 2D matrix
  X_shape <- two.d.array(raw_shape_array)
  
  # 2. Extract the model's multi-dimensional beta coefficients matrix
  all_coefs <- model_obj$coefficients
  
  # 3. Locate the index row matching our continuous 'log_size' predictor
  log_size_row <- grep("log_size", rownames(all_coefs))
  
  # Safety backup: if matrix column tracking names differ, default to index 5
  if(length(log_size_row) == 0) log_size_row <- 5 
  
  beta_slope <- all_coefs[log_size_row, ]
  
  # 4. Normalize the extracted slope vector to a unit length of 1
  beta_unit_vector <- beta_slope / sqrt(sum(beta_slope^2))
  
  # 5. Execute matrix dot product multiplication to project cases onto the vector axis
  reg_scores <- X_shape %*% beta_unit_vector
  return(as.numeric(reg_scores))
}

# Execute clean, silent matrix math calculations on both anatomical views
ant_clean_scores <- compute_silent_reg_score(ant_allom_model, ant_sym_shape)
lat_clean_scores <- compute_silent_reg_score(lat_allom_model, lat_shape)

# Assemble an unconfounded data table for ggplot2 mapping
allom_plot_df <- tibble(
  Ant_RegScore = ant_clean_scores,
  Lat_RegScore = lat_clean_scores,
  Ant_LogSize  = log(ant_csize),
  Lat_LogSize  = log(lat_csize),
  Group        = factor(ifelse(master$CPD == 1, "CPD (High Risk)", "Control (Vaginal)"))
)
# Update column designators for the plotting data frame
allom_plot_df_clean <- allom_plot_df %>%
  mutate(Group = factor(Group, levels = c("Control (Vaginal)", "CPD (High Risk)")))


# Standardize titles and themes for panel grid alignment
gg_allom_ant <- ggplot(allom_plot_df_clean, aes(x = Ant_LogSize, y = Ant_RegScore, color = Group)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("deepskyblue4", "darkorange")) +
  labs(
    title = "A. Anterior Facial Profile",
    x = "Facial Scale: log(Centroid Size)",
    y = "Shape Regression Score"
  ) +
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 11))

gg_allom_lat <- ggplot(allom_plot_df_clean, aes(x = Lat_LogSize, y = Lat_RegScore, color = Group)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("deepskyblue4", "darkorange")) +
  labs(
    title = "B. Lateral Neurocranial Profile",
    x = "Vault Scale: log(Centroid Size)",
    y = NULL
  ) +
  theme_minimal() + 
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(face = "bold", size = 11))






