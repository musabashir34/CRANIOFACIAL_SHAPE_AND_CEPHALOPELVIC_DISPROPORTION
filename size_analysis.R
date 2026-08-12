# size_analysis.R

# 1. Reconstruct and Validate the Combined Data Tables Natively
# This mirrors the clean snake_case layouts from your MANCOVA block
ancova_df_ant <- master %>%
  mutate(
    Centroid_Size   = as.numeric(ant_csize), # Ingests univariate scale vectors
    Maternal_Age    = as.numeric(`Maternal age`),
    Maternal_Height = as.numeric(`M Height(cm)`),
    Parity          = as.numeric(Parity),
    Sex             = factor(SEX, levels = c("Female", "Male")),
    Group           = factor(group, levels = c("C", "S")) # Clinical factor
  )

ancova_df_lat <- master %>%
  mutate(
    Centroid_Size   = as.numeric(lat_csize), # Ingests univariate scale vectors
    Maternal_Age    = as.numeric(`Maternal age`),
    Maternal_Height = as.numeric(`M Height(cm)`),
    Parity          = as.numeric(Parity),
    Sex             = factor(SEX, levels = c("Female", "Male")),
    Group           = factor(group, levels = c("C", "S")) # Clinical factor
  )

## ===================== A. ANTERIOR VIEW SIZE ANCOVA =====================

# Fit the sequential linear regression model (Type I Sequential Sum of Squares)
ant_size_model <- lm(
  Centroid_Size ~ Maternal_Age + Maternal_Height + Parity + Sex + Group,
  data = ancova_df_ant
)

# Extract sequential Type I ANOVA parameter frames
ant_ancova_summary <- anova(ant_size_model)

## ===================== B. LATERAL VIEW SIZE ANCOVA ======================

# Fit the sequential linear regression model (Type I Sequential Sum of Squares)
lat_size_model <- lm(
  Centroid_Size ~ Maternal_Age + Maternal_Height + Parity + Sex + Group,
  data = ancova_df_lat
)

# Extract sequential Type I ANOVA parameter frames
lat_ancova_summary <- anova(lat_size_model)

# Cache final statistics for dynamic R Markdown inline text integration
ant_size_df1 <- ant_ancova_summary["Group", "Df"]
ant_size_df2 <- ant_ancova_summary["Residuals", "Df"]
ant_size_f   <- ant_ancova_summary["Group", "F value"]
ant_size_p   <- ant_ancova_summary["Group", "Pr(>F)"]

lat_size_df1 <- lat_ancova_summary["Group", "Df"]
lat_size_df2 <- lat_ancova_summary["Residuals", "Df"]
lat_size_f   <- lat_ancova_summary["Group", "F value"]
lat_size_p   <- lat_ancova_summary["Group", "Pr(>F)"]

#===============================================================
# DISPLAY SIZE ANCOVA AS A TABLE
#===============================================================
# 1. Extract Summary Matrices Natively from the linear models
ant_summary_df <- as.data.frame(anova(ant_size_model))
lat_summary_df <- as.data.frame(anova(lat_size_model))

# 2. Re-map rows to align parameters perfectly into column slices
# Predictor rows are identical due to your matching sequential design formula
predictor_rows <- c("Maternal Age", "Maternal Height", "Parity", "Infant Sex", "CPD Group Grouping", "Residuals")

# Helper function to format p-values cleanly for journal layouts
format_p_value <- function(p_val) {
  if (is.na(p_val)) return("—")
  if (p_val < 0.001) return("< 0.001*")
  if (p_val < 0.05) return(sprintf("%.4f*", p_val))
  return(sprintf("%.4f", p_val))
}

format_f_value <- function(f_val) {
  if (is.na(f_val)) return("—")
  return(sprintf("%.2f", f_val))
}

# 3. Assemble the Side-by-Side Presentation Table
side_by_side_ancova_table <- tibble(
  `Predictor Source` = predictor_rows,
  
  # Anterior View Data Columns
  `Ant Df`    = ant_summary_df$Df,
  `Ant Sum Sq` = round(ant_summary_df$`Sum Sq`, 4),
  `Ant F-value` = sapply(ant_summary_df$`F value`, format_f_value),
  `Ant p-value` = sapply(ant_summary_df$`Pr(>F)`, format_p_value),
  
  # Lateral View Data Columns
  `Lat Df`    = lat_summary_df$Df,
  `Lat Sum Sq` = round(lat_summary_df$`Sum Sq`, 4),
  `Lat F-value` = sapply(lat_summary_df$`F value`, format_f_value),
  `Lat p-value` = sapply(lat_summary_df$`Pr(>F)`, format_p_value)
)
