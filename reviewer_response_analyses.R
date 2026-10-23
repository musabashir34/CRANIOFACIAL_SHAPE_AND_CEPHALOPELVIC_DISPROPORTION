# reviewer_response_analyses.R
# ---------------------------------------------------------------------------
# Additional analyses added in response to peer review. Source AFTER
# allometry_controlled_shape_analysis.R (needs ant_sym_shape, lat_shape, ant_csize, lat_csize, ant_group,
# lat_group, sex_numeric, N_PERM, master, ant_pca, lat_pca, mancova_df_ant, mancova_df_lat).
#
#  B. group x log(centroid size) interaction (test of allometric slope homogeneity)
#  C. sensitivity of the size-adjusted CPD group effect to the choice of size covariate
#     (log centroid size vs directly measured HC, OFD, BPD)
#  E. principal-component retention: variance explained and sensitivity of the MANOVA group effect
#
# ---------------------------------------------------------------------------

# ---- B. Interaction test -------------------------------------------------------------------------
ant_gdf_interact <- geomorph.data.frame(
  shape = ant_sym_shape, group = factor(ant_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, log_size = log(ant_csize))
lat_gdf_interact <- geomorph.data.frame(
  shape = lat_shape, group = factor(lat_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, log_size = log(lat_csize))

ant_interaction_model <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + log_size * group,
                                  data = ant_gdf_interact, iter = N_PERM, print.progress = FALSE)
lat_interaction_model <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + log_size * group,
                                  data = lat_gdf_interact, iter = N_PERM, print.progress = FALSE)
ant_int_tab <- as.data.frame(anova(ant_interaction_model)$table)
lat_int_tab <- as.data.frame(anova(lat_interaction_model)$table)
interaction_results <- data.frame(
  view = c("Anterior", "Lateral"),
  F    = c(ant_int_tab["log_size:group", "F"], lat_int_tab["log_size:group", "F"]),
  p    = c(ant_int_tab["log_size:group", "Pr(>F)"], lat_int_tab["log_size:group", "Pr(>F)"]))
write.csv(interaction_results, "interaction_group_x_logsize.csv", row.names = FALSE)

# ---- C. Sensitivity to the size covariate --------------------------------------------------------
sensitivity_results <- data.frame(size_covariate = character(), view = character(),
                                  F = numeric(), p = numeric(), stringsAsFactors = FALSE)

# --- log(Centroid Size) : Anterior ---
gdf_s <- geomorph.data.frame(shape = ant_sym_shape, group = factor(ant_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, size_cov = log(ant_csize))
m_s <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + size_cov + group,
                data = gdf_s, iter = N_PERM, print.progress = FALSE)
tab_s <- as.data.frame(anova(m_s)$table)
sensitivity_results <- rbind(sensitivity_results,
  data.frame(size_covariate = "log(Centroid Size)", view = "Anterior", F = tab_s["group", "F"], p = tab_s["group", "Pr(>F)"]))

# --- log(Centroid Size) : Lateral ---
gdf_s <- geomorph.data.frame(shape = lat_shape, group = factor(lat_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, size_cov = log(lat_csize))
m_s <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + size_cov + group,
                data = gdf_s, iter = N_PERM, print.progress = FALSE)
tab_s <- as.data.frame(anova(m_s)$table)
sensitivity_results <- rbind(sensitivity_results,
  data.frame(size_covariate = "log(Centroid Size)", view = "Lateral", F = tab_s["group", "F"], p = tab_s["group", "Pr(>F)"]))

# --- Head Circumference : Anterior ---
gdf_s <- geomorph.data.frame(shape = ant_sym_shape, group = factor(ant_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, size_cov = as.numeric(master$`HC(cm)`))
m_s <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + size_cov + group,
                data = gdf_s, iter = N_PERM, print.progress = FALSE)
tab_s <- as.data.frame(anova(m_s)$table)
sensitivity_results <- rbind(sensitivity_results,
  data.frame(size_covariate = "Head Circumference", view = "Anterior", F = tab_s["group", "F"], p = tab_s["group", "Pr(>F)"]))

# --- Head Circumference : Lateral ---
gdf_s <- geomorph.data.frame(shape = lat_shape, group = factor(lat_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, size_cov = as.numeric(master$`HC(cm)`))
m_s <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + size_cov + group,
                data = gdf_s, iter = N_PERM, print.progress = FALSE)
tab_s <- as.data.frame(anova(m_s)$table)
sensitivity_results <- rbind(sensitivity_results,
  data.frame(size_covariate = "Head Circumference", view = "Lateral", F = tab_s["group", "F"], p = tab_s["group", "Pr(>F)"]))

# --- OFD : Anterior ---
gdf_s <- geomorph.data.frame(shape = ant_sym_shape, group = factor(ant_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, size_cov = as.numeric(master$`OFD(mm)`))
m_s <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + size_cov + group,
                data = gdf_s, iter = N_PERM, print.progress = FALSE)
tab_s <- as.data.frame(anova(m_s)$table)
sensitivity_results <- rbind(sensitivity_results,
  data.frame(size_covariate = "OFD", view = "Anterior", F = tab_s["group", "F"], p = tab_s["group", "Pr(>F)"]))

# --- OFD : Lateral ---
gdf_s <- geomorph.data.frame(shape = lat_shape, group = factor(lat_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, size_cov = as.numeric(master$`OFD(mm)`))
m_s <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + size_cov + group,
                data = gdf_s, iter = N_PERM, print.progress = FALSE)
tab_s <- as.data.frame(anova(m_s)$table)
sensitivity_results <- rbind(sensitivity_results,
  data.frame(size_covariate = "OFD", view = "Lateral", F = tab_s["group", "F"], p = tab_s["group", "Pr(>F)"]))

# --- BPD : Anterior ---
gdf_s <- geomorph.data.frame(shape = ant_sym_shape, group = factor(ant_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, size_cov = as.numeric(master$`BPD(mm)`))
m_s <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + size_cov + group,
                data = gdf_s, iter = N_PERM, print.progress = FALSE)
tab_s <- as.data.frame(anova(m_s)$table)
sensitivity_results <- rbind(sensitivity_results,
  data.frame(size_covariate = "BPD", view = "Anterior", F = tab_s["group", "F"], p = tab_s["group", "Pr(>F)"]))

# --- BPD : Lateral ---
gdf_s <- geomorph.data.frame(shape = lat_shape, group = factor(lat_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, size_cov = as.numeric(master$`BPD(mm)`))
m_s <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + size_cov + group,
                data = gdf_s, iter = N_PERM, print.progress = FALSE)
tab_s <- as.data.frame(anova(m_s)$table)
sensitivity_results <- rbind(sensitivity_results,
  data.frame(size_covariate = "BPD", view = "Lateral", F = tab_s["group", "F"], p = tab_s["group", "Pr(>F)"]))

write.csv(sensitivity_results, "sensitivity_alternative_size_covariates.csv", row.names = FALSE)

# ---- E. PC retention -----------------------------------------------------------------------------
# gm.prcomp objects expose $sdev and $x (there is no $importance element in this geomorph version).
ant_cum_var <- cumsum(ant_pca$sdev^2 / sum(ant_pca$sdev^2))
lat_cum_var <- cumsum(lat_pca$sdev^2 / sum(lat_pca$sdev^2))
pc_var_summary <- data.frame(
  view = c("Anterior", "Anterior", "Lateral", "Lateral"),
  n_pc = c(5, 10, 10, 16),
  cumulative_variance = c(ant_cum_var[5], ant_cum_var[10], lat_cum_var[10], lat_cum_var[16]))

rerun_mancova_pc_count <- function(pca_obj, n_pc, mancova_df, label) {
  scores <- pca_obj$x[, 1:n_pc]
  fit <- manova(scores ~ Maternal_Age + Maternal_Height + Parity + Sex + Group, data = mancova_df)
  s <- summary(fit, test = "Pillai")
  data.frame(view = label, n_pc = n_pc, pillai = s$stats["Group", "Pillai"],
             F = s$stats["Group", "approx F"], p = s$stats["Group", "Pr(>F)"])
}
pc_sensitivity <- rbind(
  rerun_mancova_pc_count(ant_pca, 5,  mancova_df_ant, "Anterior, 5 PC (saturation-based)"),
  rerun_mancova_pc_count(ant_pca, 10, mancova_df_ant, "Anterior, 10 PC (primary)"),
  rerun_mancova_pc_count(lat_pca, 16, mancova_df_lat, "Lateral, 16 PC (saturation-based)"),
  rerun_mancova_pc_count(lat_pca, 10, mancova_df_lat, "Lateral, 10 PC (primary)"))
write.csv(pc_var_summary, "pc_variance_explained.csv", row.names = FALSE)
write.csv(pc_sensitivity, "pc_retention_sensitivity.csv", row.names = FALSE)

# ---- F. Confidence intervals on the allometric slopes --------------
# Uses allom_plot_df_clean (built in allometry_controlled_shape_analysis.R) -- the same per-specimen
# RegScore vs log(Centroid Size) data frame that draws Figure 3. The "slope" is the slope of this
# univariate shape-regression score on log(size), i.e. exactly the line geom_smooth(method="lm")
# already draws in Figure 3; this section just attaches a number and a CI to that line, per group
# and for the common within-group slope (consistent with the non-significant interaction in Section B).
# FIXED (revision 2): the common slope is now estimated with a group term (RegScore ~ LogSize + Group),
# i.e. parallel lines with separate intercepts. The earlier "pooled" row regressed RegScore on LogSize
# ignoring group, which mixes the between-group offset into the slope (anterior 0.12, below both group
# slopes). coef(fit)[2] is still the LogSize coefficient because LogSize is entered first.

slope_ci_row <- function(fit, label) {
  cf <- confint(fit)
  data.frame(model = label,
             slope = coef(fit)[2],
             ci_low = cf[2, 1],
             ci_high = cf[2, 2])
}

df_ctrl <- allom_plot_df_clean[allom_plot_df_clean$Group == "Control (Vaginal)", ]
df_case <- allom_plot_df_clean[allom_plot_df_clean$Group == "CPD cases", ]

allometric_slope_cis <- rbind(
  slope_ci_row(lm(Ant_RegScore ~ Ant_LogSize, data = df_ctrl),              "Anterior, Control"),
  slope_ci_row(lm(Ant_RegScore ~ Ant_LogSize, data = df_case),              "Anterior, CPD cases"),
  slope_ci_row(lm(Ant_RegScore ~ Ant_LogSize + Group, data = allom_plot_df_clean),  "Anterior, common within-group slope"),
  slope_ci_row(lm(Lat_RegScore ~ Lat_LogSize, data = df_ctrl),              "Lateral, Control"),
  slope_ci_row(lm(Lat_RegScore ~ Lat_LogSize, data = df_case),              "Lateral, CPD cases"),
  slope_ci_row(lm(Lat_RegScore ~ Lat_LogSize + Group, data = allom_plot_df_clean),  "Lateral, common within-group slope")
)
rownames(allometric_slope_cis) <- NULL

cat("\n=== Allometric slope estimates with 95% CIs ===\n")
print(allometric_slope_cis)
write.csv(allometric_slope_cis, "allometric_slope_cis.csv", row.names = FALSE)

# ---- G. Gestational age as a covariate --------------------------------
# GA_weeks is merged into master by data_ingestion.R (completed weeks; all neonates term, 38-41 wk).
# (i) Logistic regression: Table 1 model with GA added (reuses cpd_glm / reg_df from
#     traditional_anthropometric_regression.R).
reg_df$GA_weeks <- as.numeric(master$GA_weeks)
cpd_glm_ga <- update(cpd_glm, . ~ . + GA_weeks)
co_ga <- summary(cpd_glm_ga)$coefficients[-1, , drop = FALSE]
ga_regression <- data.frame(
  Variable    = rownames(co_ga),
  aOR         = round(exp(co_ga[, "Estimate"]), 2),
  Lower_95_CI = round(exp(co_ga[, "Estimate"] - z_crit * co_ga[, "Std. Error"]), 2),
  Upper_95_CI = round(exp(co_ga[, "Estimate"] + z_crit * co_ga[, "Std. Error"]), 2),
  p_value     = round(co_ga[, "Pr(>|z|)"], 4), row.names = NULL)
ga_lrt <- anova(cpd_glm, cpd_glm_ga, test = "Chisq")
cat("\n=== Logistic regression with gestational age (per week) added ===\n"); print(ga_regression)
cat("LR test for adding GA: chi2 =", round(ga_lrt$Deviance[2], 3), ", p =", round(ga_lrt$`Pr(>Chi)`[2], 3), "\n")
write.csv(ga_regression, "table_S6_logistic_with_gestational_age.csv", row.names = FALSE)

# (ii) Shape models: GA entered with the other covariates (before log size and group, as in the main model)
ant_gdf_ga <- geomorph.data.frame(
  shape = ant_sym_shape, group = factor(ant_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, GA = as.numeric(master$GA_weeks),
  log_size = log(ant_csize))
lat_gdf_ga <- geomorph.data.frame(
  shape = lat_shape, group = factor(lat_group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, GA = as.numeric(master$GA_weeks),
  log_size = log(lat_csize))
ant_ga_model <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + GA + log_size + group,
                         data = ant_gdf_ga, iter = N_PERM, print.progress = FALSE)
lat_ga_model <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + GA + log_size + group,
                         data = lat_gdf_ga, iter = N_PERM, print.progress = FALSE)
ant_ga_tab <- as.data.frame(anova(ant_ga_model)$table)
lat_ga_tab <- as.data.frame(anova(lat_ga_model)$table)
ga_shape_results <- data.frame(
  view = c("Anterior", "Lateral"),
  GA_F = c(ant_ga_tab["GA", "F"], lat_ga_tab["GA", "F"]),
  GA_p = c(ant_ga_tab["GA", "Pr(>F)"], lat_ga_tab["GA", "Pr(>F)"]),
  group_F = c(ant_ga_tab["group", "F"], lat_ga_tab["group", "F"]),
  group_p = c(ant_ga_tab["group", "Pr(>F)"], lat_ga_tab["group", "Pr(>F)"]))
cat("\n=== Shape models with gestational age added: CPD group effect ===\n"); print(ga_shape_results)
cat("Original (no GA): anterior F = 5.36, p < 0.001; lateral F = 3.83, p = 0.0075\n")
write.csv(ga_shape_results, "gestational_age_shape_models.csv", row.names = FALSE)
