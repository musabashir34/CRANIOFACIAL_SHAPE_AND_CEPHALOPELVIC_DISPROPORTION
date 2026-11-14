# check_lateral_scale.R
# ---------------------------------------------------------------------------
# Diagnostic for the two clusters in lateral log(centroid size) (Figure 3B:
# ~7.25 and ~8.05, a ~2.2-fold difference in pixel size).
# Source AFTER the main pipeline (needs lat_land, lat_codes, lat_csize, lat_shape,
# master, sex_numeric, lat_group, N_PERM).
# ---------------------------------------------------------------------------

lat_logcs <- log(lat_csize)
cut_pt    <- 7.6                       # gap between the two clusters in Figure 3B
lat_cluster <- factor(ifelse(lat_logcs > cut_pt, "large", "small"))

cat("\n== Cluster sizes by group ==\n")
print(table(cluster = lat_cluster, group = lat_codes$group))

cat("\n== Specimens in each cluster (look for a pattern in IDs / dates / file names) ==\n")
print(split(dimnames(lat_land)[[3]], lat_cluster))

# 1. Do some TPS records carry a SCALE= line (mm) and others not (pixels)?
tps_lines <- readLines("LATERAL_IMAGES_4_SAMPLE_AND_CONTROL.TPS")
img_idx   <- grep("^IMAGE=", tps_lines)
scale_at  <- sapply(seq_along(img_idx), function(i) {
  start <- if (i == 1) 1 else img_idx[i - 1] + 1
  s <- grep("^SCALE=", tps_lines[start:img_idx[i]], value = TRUE)
  if (length(s)) sub("SCALE=", "", s[1]) else NA
})
cat("\n== SCALE= lines in the lateral TPS file ==\n")
print(table(has_scale = !is.na(scale_at), cluster = lat_cluster, useNA = "ifany"))

# 2. Do the image files differ in pixel dimensions? (needs the 'jpeg' package and the image files)
 img_files <- sub("IMAGE=", "", tps_lines[img_idx])
 dims <- t(sapply(img_files, function(f) dim(jpeg::readJPEG(f))[1:2]))
 print(table(paste(dims[, 1], "x", dims[, 2]), lat_cluster))

# 3. Does the lateral group effect survive with cluster as a covariate?
lat_gdf_cl <- geomorph.data.frame(
  shape = lat_shape, group = factor(lat_group), cluster = lat_cluster,
  maternal_age = as.numeric(master$`Maternal age`), maternal_ht = as.numeric(master$`M Height(cm)`),
  parity = as.numeric(master$Parity), sex = sex_numeric, log_size = lat_logcs)
lat_cl_model <- procD.lm(shape ~ maternal_age + maternal_ht + parity + sex + cluster + log_size + group,
                         data = lat_gdf_cl, iter = N_PERM, print.progress = FALSE)
cat("\n== Lateral shape model with acquisition cluster as covariate ==\n")
print(anova(lat_cl_model)$table)

# 4. Lateral size ANCOVA with cluster as covariate
lat_size_cl <- lm(lat_csize ~ as.numeric(master$`Maternal age`) + as.numeric(master$`M Height(cm)`) +
                    as.numeric(master$Parity) + factor(master$SEX) + lat_cluster + factor(lat_group))
cat("\n== Lateral centroid-size ANCOVA with cluster ==\n")
print(anova(lat_size_cl))
