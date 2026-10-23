# traditional_anthropometric_regression.R
# ---------------------------------------------------------------------------
# Multivariable logistic regression of CPD case status on neonatal anthropometry and maternal covariates
# (Table 1). Requires the data frame `master` built by build_master_dataset_geomorph.R.
#
# NOTE: this file REPLACES the earlier version of the script. Changes:
#   * Birth weight is modelled in kilograms (SI). The raw spreadsheet column `Weight(ibs)` is in POUNDS
#     (1 lb = 0.45359237 kg); the earlier Table 1 reported the per-POUND odds ratio (2.76) but labelled it "kg".
#   * 95% CIs are Wald intervals (estimate +/- 1.96 SE on the log-odds scale), consistent with the Wald p-values.
#   * Variance inflation factors are reported (computed from the coefficient covariance matrix, i.e. the
#     definition used by car::vif() for GLMs; no extra package required).
#   * Odds ratios are per 1 unit of each predictor (unit stated in the row label used in the manuscript).
# Only base R (stats) is used.
# ---------------------------------------------------------------------------

lb_to_kg <- 0.45359237

reg_df <- data.frame(
  CPD             = as.numeric(master$CPD),
  OFD_mm          = as.numeric(master$`OFD(mm)`),
  Weight_kg       = as.numeric(master$`Weight(ibs)`) * lb_to_kg,
  Maternal_Age    = as.numeric(master$`Maternal age`),
  HC_cm           = as.numeric(master$`HC(cm)`),
  BPD_mm          = as.numeric(master$`BPD(mm)`),
  Maternal_Height = as.numeric(master$`M Height(cm)`),
  Parity          = as.numeric(master$Parity),
  SexMale         = as.numeric(trimws(as.character(master$SEX)) == "Male")
)

# Simultaneous entry of all predictors (no stepwise selection)
cpd_glm <- glm(CPD ~ OFD_mm + Weight_kg + Maternal_Age + HC_cm + BPD_mm +
                 Maternal_Height + Parity + SexMale,
               data = reg_df, family = binomial(link = "logit"))

coef_tab <- summary(cpd_glm)$coefficients
coef_tab <- coef_tab[rownames(coef_tab) != "(Intercept)", , drop = FALSE]
z_crit   <- qnorm(0.975)

clinical_table <- data.frame(
  Variable     = rownames(coef_tab),
  aOR          = round(exp(coef_tab[, "Estimate"]), 2),
  Lower_95_CI  = round(exp(coef_tab[, "Estimate"] - z_crit * coef_tab[, "Std. Error"]), 2),
  Upper_95_CI  = round(exp(coef_tab[, "Estimate"] + z_crit * coef_tab[, "Std. Error"]), 2),
  p_value      = round(coef_tab[, "Pr(>|z|)"], 4),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Variance inflation factors (car::vif definition for GLMs: from the coefficient covariance matrix)
vif_from_vcov <- function(fit) {
  v <- vcov(fit)
  v <- v[rownames(v) != "(Intercept)", colnames(v) != "(Intercept)", drop = FALSE]
  r <- cov2cor(v)
  setNames(diag(solve(r)), rownames(r))
}
vif_values <- vif_from_vcov(cpd_glm)
clinical_table$VIF <- round(vif_values[clinical_table$Variable], 2)

# Model fit
ll_full <- as.numeric(logLik(cpd_glm))
ll_null <- as.numeric(logLik(update(cpd_glm, . ~ 1)))
n_obs   <- nrow(reg_df)
r2_cox_snell <- 1 - exp((2 / n_obs) * (ll_null - ll_full))
r2_nagelkerke <- r2_cox_snell / (1 - exp((2 / n_obs) * ll_null))
lr_chisq <- cpd_glm$null.deviance - cpd_glm$deviance
lr_df    <- cpd_glm$df.null - cpd_glm$df.residual
lr_p     <- pchisq(lr_chisq, lr_df, lower.tail = FALSE)

# Leave-one-out cross-validated classification accuracy (descriptive discrimination in this
# artificially balanced 1:1 case-control sample -- NOT a measure of prospective screening performance)
loocv_pred <- rep(NA_real_, n_obs)
for (i in seq_len(n_obs)) {
  fit_i <- suppressWarnings(glm(formula(cpd_glm), data = reg_df[-i, ], family = binomial))
  loocv_pred[i] <- predict(fit_i, newdata = reg_df[i, , drop = FALSE], type = "response")
}
loocv_accuracy <- mean((loocv_pred >= 0.5) == (reg_df$CPD == 1))

# Convenience objects used in the manuscript text
pick <- function(v) clinical_table[clinical_table$Variable == v, ]
OFD_OR <- pick("OFD_mm")$aOR; OFD_LCI <- pick("OFD_mm")$Lower_95_CI; OFD_UCI <- pick("OFD_mm")$Upper_95_CI; OFD_p <- pick("OFD_mm")$p_value
BW_OR  <- pick("Weight_kg")$aOR; BW_LCI <- pick("Weight_kg")$Lower_95_CI; BW_UCI <- pick("Weight_kg")$Upper_95_CI; BW_p <- pick("Weight_kg")$p_value
BPD_p  <- pick("BPD_mm")$p_value
HC_p   <- pick("HC_cm")$p_value

# Birth weight per 100 g (clinically interpretable increment), Wald CI
bw_b  <- coef_tab["Weight_kg", "Estimate"]; bw_se <- coef_tab["Weight_kg", "Std. Error"]
BW_100g_OR  <- round(exp(bw_b / 10), 2)
BW_100g_LCI <- round(exp((bw_b - z_crit * bw_se) / 10), 2)
BW_100g_UCI <- round(exp((bw_b + z_crit * bw_se) / 10), 2)

write.csv(clinical_table, "table1_logistic_regression.csv", row.names = FALSE)
