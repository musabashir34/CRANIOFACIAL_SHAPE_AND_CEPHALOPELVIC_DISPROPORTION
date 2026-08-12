################################################################################
# traditional_anthropometric_regression.R
#
# Runs simultaneous multivariable binary logistic regression on traditional
# anthropometric variables and maternal confounders 
# REQUIRES: dplyr, broom, car, pROC
# install.packages(c("dplyr", "broom", "car", "pROC"))
################################################################################

# 1. Structure the dataset for medical regression logic
regression_data <- master %>%
  mutate(
    # Ensure your target factor is binary (1 = Emergency C-Section for CPD, 0 = Vaginal Control)
    CPD = factor(CPD, levels = c(0,1)),
    Maternal_Age = as.numeric(master$`Maternal age`),
    Maternal_Height = as.numeric(master$`M Height(cm)`),
    Parity = as.numeric(Parity),
    Sex = factor(SEX, levels = c("Female", "Male"))
  )

# 2. Fit the Simultaneous Binary Logistic Regression Model (Type III Logic)
# Entering all predictors into a single formula forces R to execute a 
# maximum-likelihood simultaneous estimation. Every predictor is evaluated 
# while holding all other infant and maternal variables completely constant.
anthropometric_model <- glm(
  CPD ~ `OFD(mm)` + 
    `Weight(ibs)` + 
    Maternal_Age + 
    `HC(cm)` + 
    `BPD(mm)` + 
    Maternal_Height + 
    Parity +
    Sex, 
  data = regression_data,
  family = binomial(link = "logit")
)

# 3. Extract and format the Adjusted Odds Ratios (aOR) & 95% Confidence Intervals
model_summary <- summary(anthropometric_model)

clinical_table <- tidy(anthropometric_model, exponentiate = FALSE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    # Exponentiate the raw log-odds estimates to calculate Adjusted Odds Ratios (aOR)
    aOR = exp(estimate),
    # Compute robust Profile Likelihood Confidence Intervals
    Lower_95_CI = exp(estimate - (1.96 * std.error)),
    Upper_95_CI = exp(estimate + (1.96 * std.error)),
    # Format p-values for standard text rendering
    p_value = round(p.value, 4)
  ) %>%
  select(Variable = term, aOR, Lower_95_CI, Upper_95_CI, p_value)


# 4. Diagnostic Model Fit Statistics
# A. Nagelkerke's Pseudo-R2 (Variance explanation estimation)
null_model <- glm(CPD ~ 1, data = regression_data, family = binomial(link = "logit"))
log_L_model <- as.numeric(logLik(anthropometric_model))
log_L_null  <- as.numeric(logLik(null_model))
n_samples   <- nrow(regression_data)

r2_cox_snell <- 1 - exp((2 / n_samples) * (log_L_null - log_L_model))
r2_nagelkerke <- r2_cox_snell / (1 - exp((2 / n_samples) * log_L_null))

# B. Collinearity Diagnostic: Variance Inflation Factor (VIF)
# Anthropometric variables are naturally highly correlated. VIF checks if severe 
# multicollinearity is destabilizing the regression weight estimations.
vif_values <- car::vif(anthropometric_model)

# Extract values for inline r code
OFD_OR <- round(clinical_table$aOR[clinical_table$Variable == "`OFD(mm)`"], 2)
OFD_LCI <- round(clinical_table$Lower_95_CI[clinical_table$Variable == "`OFD(mm)`"], 2)
OFD_UCI <- round(clinical_table$Upper_95_CI[clinical_table$Variable == "`OFD(mm)`"], 2)
OFD_p <- clinical_table$p_value[clinical_table$Variable == "`OFD(mm)`"]
BW_OR <- round(clinical_table$aOR[clinical_table$Variable == "`Weight(ibs)`"], 2)
BW_LCI <- round(clinical_table$Lower_95_CI[clinical_table$Variable == "`Weight(ibs)`"], 2)
BW_UCI <- round(clinical_table$Upper_95_CI[clinical_table$Variable == "`Weight(ibs)`"], 2)
BW_p <- clinical_table$p_value[clinical_table$Variable == "`Weight(ibs)`"]
BPD_p <- clinical_table$p_value[clinical_table$Variable == "`BPD(mm)`"]
HC_p <- clinical_table$p_value[clinical_table$Variable == "`HC(cm)`"]
