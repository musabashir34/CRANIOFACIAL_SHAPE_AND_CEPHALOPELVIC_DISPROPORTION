## ============================================================================
## build_master_dataset_geomorph.R  (geomorph/Morpho version)
##
## Builds the covariate + GM master data frame from the geomorph-derived PCA
## scores and centroid sizes
## ============================================================================


N_ANT_PC <- 5
N_LAT_PC <- 16

ant_shape_df <- data.frame(
  group = ant_codes$group, id = ant_codes$id,
  A_CS = ant_csize
)
for (i in 1:N_ANT_PC) ant_shape_df[[paste0("A_symPC", i)]] <- ant_pca$x[, i]

lat_shape_df <- data.frame(
  group = lat_codes$group, id = lat_codes$id,
  L_CS = lat_csize
)
for (i in 1:N_LAT_PC) lat_shape_df[[paste0("L_PC", i)]] <- lat_pca$x[, i]

master <- covar %>%
  left_join(ant_shape_df, by = c("group", "id")) %>%
  left_join(lat_shape_df, by = c("group", "id"))

names(master) <- trimws(names(master))

write.csv(master, "master_linked_dataset_geomorph.csv", row.names = FALSE)
