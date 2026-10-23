## ============================================================================
## landmark_exclusion_sensitivity.R
##
## This is a STANDALONE script. It does not modify data_ingestion.R,
## primary_gm_analysis_geomorph.R, or any other existing script, and it does
## not overwrite ant_sym_shape / lat_shape / ant_csize / lat_csize / ant_pca /
## lat_pca or any other object the main results depend on -- every object
## created here has a _excl suffix. This keeps the sensitivity check fully
## separate from, and incapable of silently altering, the primary analysis.
##
## RUN ORDER: source this AFTER data_ingestion.R, primary_gm_analysis_geomorph.R,
## and build_master_dataset_geomorph.R (needs ant_land, lat_land, ant_ind,
## ant_codes, lat_codes, master, N_PERM, sex_numeric).
##
## Landmark indices (from the Figure 1 / Figure 2 legends and the land.pairs
## matrix in primary_gm_analysis_geomorph.R):
##   Anterior: zygion = landmarks 14 (left) and 15 (right) -- conveniently the
##             LAST bilateral pair in ant_land.pairs, so dropping them does not
##             renumber any other landmark.
##   Lateral:  opisthocranion = landmark 11 (unpaired; the lateral view has no
##             bilateral symmetry step).
## ============================================================================

library(geomorph)

## ---- 1. Anterior: drop zygion (landmarks 14, 15) -----------------------------------------------

ant_land_excl <- ant_land[1:13, , ]                       # drop the trailing pair cleanly
ant_pairs_excl <- matrix(c(6,7, 8,9, 10,11, 12,13), ncol = 2, byrow = TRUE)  # same pairs, minus 14-15

ant_gpa_plain_excl <- gpagen(ant_land_excl, print.progress = FALSE, curves = NULL)
ant_csize_excl <- ant_gpa_plain_excl$Csize

bs_excl <- bilat.symmetry(A = ant_land_excl, ind = ant_ind,
                          object.sym = TRUE, land.pairs = ant_pairs_excl,
                          iter = N_PERM, print.progress = FALSE)
ant_sym_shape_excl <- bs_excl$symm.shape

## ---- 2. Lateral: drop opisthocranion (landmark 11) ---------------------------------------------

lat_land_excl <- lat_land[-11, , ]
lat_gpa_excl <- gpagen(lat_land_excl, print.progress = FALSE, curves = NULL)
lat_csize_excl <- lat_gpa_excl$Csize
lat_shape_excl <- lat_gpa_excl$coords

## ---- 3. Re-run the confounder + allometry group-effect test on the reduced landmark sets --------
## Same model formula and covariate order as allometry_controlled_shape_analysis.R
## (shape ~ maternal_age + maternal_ht + parity + sex + log_size + group), so the CPD-group row is
## directly comparable to the original, full-landmark result (anterior F = 5.36, p < 0.001; lateral
## F = 3.83, p = 0.0075).

ant_gdf_excl <- geomorph.data.frame(
  shape = ant_sym_shape_excl, group = factor(ant_codes$group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, log_size = log(ant_csize_excl))

lat_gdf_excl <- geomorph.data.frame(
  shape = lat_shape_excl, group = factor(lat_codes$group),
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, log_size = log(lat_csize_excl))

ant_allom_model_excl <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + log_size + group,
                                 data = ant_gdf_excl, iter = N_PERM, print.progress = FALSE)
lat_allom_model_excl <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + log_size + group,
                                 data = lat_gdf_excl, iter = N_PERM, print.progress = FALSE)

ant_tab_excl <- as.data.frame(anova(ant_allom_model_excl)$table)
lat_tab_excl <- as.data.frame(anova(lat_allom_model_excl)$table)

landmark_exclusion_results <- data.frame(
  view = c("Anterior (zygion excluded)", "Lateral (opisthocranion excluded)"),
  n_landmarks = c(13, 13),
  F = c(ant_tab_excl["group", "F"], lat_tab_excl["group", "F"]),
  p = c(ant_tab_excl["group", "Pr(>F)"], lat_tab_excl["group", "Pr(>F)"])
)

cat("\n=== CPD group effect with the most error-prone landmark excluded ===\n")
print(landmark_exclusion_results)
cat("\nOriginal (full landmark set), for comparison: anterior F = 5.36, p < 0.001; lateral F = 3.83, p = 0.0075.\n")

write.csv(landmark_exclusion_results, "landmark_exclusion_sensitivity.csv", row.names = FALSE)
