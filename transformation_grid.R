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

