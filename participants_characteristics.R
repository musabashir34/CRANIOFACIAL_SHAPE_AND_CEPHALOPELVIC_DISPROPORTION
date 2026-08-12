# participants_characteristics.R

cat("====================================================================\n")
cat("   GENERATING UPDATED SUPPLEMENTARY TABLE 1: ANTHROPOMETRIC SPECS  \n")
cat("====================================================================\n\n")

# 1. Clean and Prepare the Synchronized Data Matrices
supp_df <- master %>%
  mutate(
    Group           = factor(group, levels = c("C", "S"), labels = c("Control", "CPD Case")),
    Maternal_Age    = as.numeric(`Maternal age`),
    Maternal_Height = as.numeric(`M Height(cm)`),
    Parity          = as.numeric(Parity),
    Sex             = factor(SEX),
    
    # INJECTING BIOMETRIC INFANT BIOPARAMETERS (Adjust string tags if your sheet names differ)
    Birth_Weight    = as.numeric(`Weight(ibs)`),
    Head_Circ       = as.numeric(`HC(cm)`),
    Biparietal_Diam    = as.numeric(`BPD(mm)`),
    OccipitoFrontal_Diam      = as.numeric(`OFD(mm)`)
  )

# 2. Run Baseline Inferential Testing Matrices (Extracting Independent t-test p-values)
p_age      <- t.test(Maternal_Age ~ Group, data = supp_df)$p.value
p_ht       <- t.test(Maternal_Height ~ Group, data = supp_df)$p.value
p_par      <- t.test(Parity ~ Group, data = supp_df)$p.value
p_sex      <- chisq.test(table(supp_df$Sex, supp_df$Group))$p.value

p_weight   <- t.test(Birth_Weight ~ Group, data = supp_df)$p.value
p_head_c   <- t.test(Head_Circ ~ Group, data = supp_df)$p.value
p_BPD   <- t.test(Biparietal_Diam ~ Group, data = supp_df)$p.value
p_OFD  <- t.test(OccipitoFrontal_Diam ~ Group, data = supp_df)$p.value

# Helper function to format p-values for clear journal layout presentation
format_supp_p <- function(p) {
  if (is.na(p)) return("—")
  if (p < 0.001) return("< 0.001*")
  if (p < 0.05) return(sprintf("%.4f*", p))
  return(sprintf("%.4f", p))
}

# 3. Compile Means and Standard Deviations by Clinical Group
summary_stats <- supp_df %>%
  group_by(Group) %>%
  summarise(
    N_Count        = n(),
    Age_Mean       = mean(Maternal_Age, na.rm = TRUE),
    Age_SD         = sd(Maternal_Age, na.rm = TRUE),
    Height_Mean    = mean(Maternal_Height, na.rm = TRUE),
    Height_SD      = sd(Maternal_Height, na.rm = TRUE),
    Parity_Mean    = mean(Parity, na.rm = TRUE),
    Parity_SD      = sd(Parity, na.rm = TRUE),
    Male_Count     = sum(Sex == "Male", na.rm = TRUE),
    Female_Count   = sum(Sex == "Female", na.rm = TRUE),
    
    # Bioparameter Summaries
    Weight_Mean    = mean(Birth_Weight, na.rm = TRUE),
    Weight_SD      = sd(Birth_Weight, na.rm = TRUE),
    HeadC_Mean     = mean(Head_Circ, na.rm = TRUE),
    HeadC_SD       = sd(Head_Circ, na.rm = TRUE),
    BPD_Mean    = mean(Biparietal_Diam, na.rm = TRUE),
    BPD_SD      = sd(Biparietal_Diam, na.rm = TRUE),
    OFD_Mean    = mean(OccipitoFrontal_Diam, na.rm = TRUE),
    OFD_SD      = sd(OccipitoFrontal_Diam, na.rm = TRUE)
  )

# Split group frames for mapping row cell injections
ctrl <- summary_stats %>% filter(Group == "Control")
cpd  <- summary_stats %>% filter(Group == "CPD Case")

# 4. Construct the Final Styled Supplementary Frame Layout
supplementary_table_01 <- tibble(
  `Baseline Participant & Neonatal Characteristic` = c(
    "Total Sample Size (n)",
    "**Maternal Demographics**",
    "Maternal Age (Years, Mean ± SD)",
    "Maternal Height (cm, Mean ± SD)",
    "Birth Parity (Count, Mean ± SD)",
    "**Neonatal Anthropometrics**",
    "Infant Biological Sex: Male (n, %)",
    "Infant Biological Sex: Female (n, %)",
    "Birth Weight (lbs, Mean ± SD)",
    "OFC Head Circumference (cm, Mean ± SD)",
    "Biparietal Diameter (mm, Mean ± SD)",
    "Occipitofrontal Diameter (mm, Mean ± SD)"
  ),
  `Control Group (Vaginal Delivery)` = c(
    as.character(ctrl$N_Count),
    "—",
    sprintf("%.2f ± %.2f", ctrl$Age_Mean, ctrl$Age_SD),
    sprintf("%.2f ± %.2f", ctrl$Height_Mean, ctrl$Height_SD),
    sprintf("%.2f ± %.2f", ctrl$Parity_Mean, ctrl$Parity_SD),
    "—",
    sprintf("%d (%.1f%%)", ctrl$Male_Count, (ctrl$Male_Count / ctrl$N_Count) * 100),
    sprintf("%d (%.1f%%)", ctrl$Female_Count, (ctrl$Female_Count / ctrl$N_Count) * 100),
    sprintf("%.1f ± %.1f", ctrl$Weight_Mean, ctrl$Weight_SD),
    sprintf("%.2f ± %.2f", ctrl$HeadC_Mean, ctrl$HeadC_SD),
    sprintf("%.2f ± %.2f", ctrl$BPD_Mean, ctrl$BPD_SD),
    sprintf("%.2f ± %.2f", ctrl$OFD_Mean, ctrl$OFD_SD)
  ),
  `CPD Case Group (Emergency CS)` = c(
    as.character(cpd$N_Count),
    "—",
    sprintf("%.2f ± %.2f", cpd$Age_Mean, cpd$Age_SD),
    sprintf("%.2f ± %.2f", cpd$Height_Mean, cpd$Height_SD),
    sprintf("%.2f ± %.2f", cpd$Parity_Mean, cpd$Parity_SD),
    "—",
    sprintf("%d (%.1f%%)", cpd$Male_Count, (cpd$Male_Count / cpd$N_Count) * 100),
    sprintf("%d (%.1f%%)", cpd$Female_Count, (cpd$Female_Count / cpd$N_Count) * 100),
    sprintf("%.1f ± %.1f", cpd$Weight_Mean, cpd$Weight_SD),
    sprintf("%.2f ± %.2f", cpd$HeadC_Mean, cpd$HeadC_SD),
    sprintf("%.2f ± %.2f", cpd$BPD_Mean, cpd$BPD_SD),
    sprintf("%.2f ± %.2f", cpd$OFD_Mean, cpd$OFD_SD)
  ),
  `Statistical Test p-value` = c(
    "—", "—",
    format_supp_p(p_age), 
    format_supp_p(p_ht), 
    format_supp_p(p_par), 
    "—",
    format_supp_p(p_sex), 
    "Included in row above",
    format_supp_p(p_weight),
    format_supp_p(p_head_c),
    format_supp_p(p_BPD),
    format_supp_p(p_OFD)
  )
)

# 5. Export as a standalone CSV matrix file for immediate journal compilation
write.csv(supplementary_table_01, "Supplementary_Table_01_Participant_Characteristics.csv", row.names = FALSE)
cat("SUCCESS: Updated sheet file with infant anthropometrics has been exported!\n")
