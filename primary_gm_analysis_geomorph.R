## ============================================================================
## primary_gm_analysis_geomorph.R
##
## REQUIRES: geomorph, Morpho, MASS
## install.packages(c("geomorph", "Morpho"))
## ============================================================================

set.seed(42)
N_PERM <- 9999


## ================================ ANTERIOR ==================================

ant_group <- ant_codes$group
ant_ind   <- dimnames(ant_land)[[3]]

## Centroid size (computed on the raw, pre-reflection configuration --
ant_gpa_plain <- gpagen(ant_land, print.progress = FALSE, curves = NULL)
ant_csize <- ant_gpa_plain$Csize

## Object-symmetry GPA: 5 midline landmarks (1-5) + 5 bilateral pairs
## (6-7, 8-9, 10-11, 12-13, 14-15), matching the pairing validated against
ant_land.pairs <- matrix(c(6,7, 8,9, 10,11, 12,13, 14,15), ncol = 2, byrow = TRUE)

bs <- bilat.symmetry(A = ant_land, ind = ant_ind,
                     object.sym = TRUE, land.pairs = ant_land.pairs,
                     iter = N_PERM, print.progress = FALSE)

## Symmetric-component shape coordinates:
ant_sym_shape <- bs$symm.shape

ant_pca <- gm.prcomp(ant_sym_shape)

ant_gdf <- geomorph.data.frame(shape = ant_sym_shape, group = ant_group)

ant_pc_scores <- ant_pca$x[, 1:10]


## ================================= LATERAL ===================================

lat_group <- lat_codes$group
lat_gpa <- gpagen(lat_land, print.progress = FALSE, curves = NULL)
lat_csize <- lat_gpa$Csize
lat_shape <- lat_gpa$coords

lat_pca <- gm.prcomp(lat_shape)

lat_gdf <- geomorph.data.frame(shape = lat_shape, group = lat_group)
lat_pc_scores <- lat_pca$x[, 1:10]
