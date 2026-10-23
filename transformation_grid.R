# transformation_grid.R

# 1. Isolate structural landmark dimensions natively from your coordinate tracking slots
p_ant <- dim(ant_sym_shape)[1]; k_ant <- dim(ant_sym_shape)[2]
p_lat <- dim(lat_shape)[1];     k_lat <- dim(lat_shape)[2]

# 2. Extract model-fitted coordinates back into 3D landmark array structures
ant_fitted_array <- arrayspecs(ant_allom_model$fitted, p = p_ant, k = k_ant)
lat_fitted_array <- arrayspecs(lat_allom_model$fitted, p = p_lat, k = k_lat)

# 3. FIXED DATA ALIGNMENT: Sort array slices using index-matched metadata group keys
ant_C_adj <- mshape(ant_fitted_array[, , master$CPD == 0])
ant_S_adj <- mshape(ant_fitted_array[, , master$CPD == 1])

lat_C_adj <- mshape(lat_fitted_array[, , master$CPD == 0])
lat_S_adj <- mshape(lat_fitted_array[, , master$CPD == 1])


# 4. GRAPHICAL ABSTRACT EXPORT (Morphologie guide for authors)
# Minimum 531 x 1328 px (h x w), readable at 5 x 13 cm; TIFF accepted.
# 13 x 5 cm at 300 dpi = 1535 x 591 px, which meets the minimum.
if (!dir.exists("figures_publication")) dir.create("figures_publication")
tiff("figures_publication/Graphical_Abstract.tiff",
     width = 13, height = 5, units = "cm", res = 300, compression = "lzw")
layout(matrix(1:2, ncol = 2))
op_ga <- par(mar = c(0.5, 0.5, 1.8, 0.5))
plotRefToTarget(ant_C_adj, ant_S_adj, method = "TPS", mag = 3)
title(main = "Anterior view", cex.main = 0.9)
plotRefToTarget(lat_C_adj, lat_S_adj, method = "TPS", mag = 3)
title(main = "Lateral view", cex.main = 0.9)
par(op_ga)
layout(1)
dev.off()
