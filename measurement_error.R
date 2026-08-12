## ============================================================================
## measurement_error.R  (geomorph/Morpho pipeline)
## Intra-observer digitising error assessment (Klingenberg & McIntyre, 1998;
## Fruciano, 2016), using geomorph::gpagen() to superimpose the combined
## original + replicate landmark configurations for a redigitised subsample,
## then partitioning shape variance into an "Individual" term (biological
## variation) and an "Error" term (measurement error) via a two-way
## Procrustes ANOVA. Requires the replicate TPS files described in the main
## analysis report (ANTERIOR_REPLICATE.TPS, LATERAL_REPLICATE.TPS).
## ============================================================================

set.seed(42)

## ---- Build matched original + replicate arrays for the subsample ---------
build_error_input <- function(original_tps, replicate_tps) {
  orig <- readland.tps(original_tps, specID = "imageID", warnmsg = FALSE)
  repl <- readland.tps(replicate_tps, specID = "imageID", warnmsg = FALSE)
  
  parse_code <- function(names_vec) {
    m <- regmatches(names_vec, regexpr("[AL][SC][0-9]{3}", names_vec, ignore.case = TRUE, perl = TRUE))
    paste0(toupper(substr(m, 2, 2)), as.integer(substr(m, 3, 5)))
  }
  key_orig <- parse_code(dimnames(orig)[[3]])
  key_repl <- parse_code(dimnames(repl)[[3]])
  
  common <- intersect(key_repl, key_orig)
  stopifnot(length(common) == length(key_repl))
  idx_orig <- match(common, key_orig)
  idx_repl <- match(common, key_repl)
  
  k <- dim(orig)[1]; n <- length(common)
  combined <- array(NA_real_, dim = c(k, 2, 2 * n))
  combined[, , 1:n] <- orig[, , idx_orig]
  combined[, , (n + 1):(2 * n)] <- repl[, , idx_repl]
  
  list(coords = combined, individual = rep(common, 2), n_specimens = n)
}

## ---- Two-way Procrustes ANOVA: Individual vs. Error -----------------------
measurement_error_anova <- function(coords_array, individual) {
  gpa_res <- gpagen(coords_array, print.progress = FALSE, curves = NULL)
  aligned <- gpa_res$coords
  coords_2d <- two.d.array(aligned)
  N <- nrow(coords_2d)
  
  ind_factor <- factor(individual)
  n_ind <- nlevels(ind_factor)
  reps_per_ind <- N / n_ind
  
  grand_mean <- colMeans(coords_2d)
  ind_means <- t(sapply(levels(ind_factor), function(g) colMeans(coords_2d[ind_factor == g, , drop = FALSE])))
  
  ss_total <- sum(sweep(coords_2d, 2, grand_mean)^2)
  ss_individual <- reps_per_ind * sum(sweep(ind_means, 2, grand_mean)^2)
  ss_error <- ss_total - ss_individual
  
  df_individual <- n_ind - 1
  df_error <- N - n_ind
  ms_individual <- ss_individual / df_individual
  ms_error <- ss_error / df_error
  f_stat <- ms_individual / ms_error
  
  n_perm <- 9999
  perm_f <- sapply(seq_len(n_perm), function(p) {
    perm_ind <- sample(ind_factor)
    perm_means <- t(sapply(levels(perm_ind), function(g) colMeans(coords_2d[perm_ind == g, , drop = FALSE])))
    ss_ind_p <- reps_per_ind * sum(sweep(perm_means, 2, grand_mean)^2)
    ss_err_p <- ss_total - ss_ind_p
    (ss_ind_p / df_individual) / (ss_err_p / df_error)
  })
  pval <- (sum(perm_f >= f_stat) + 1) / (n_perm + 1)
  
  pct_error <- 100 * ss_error / ss_total
  pct_individual <- 100 * ss_individual / ss_total
  repeatability <- ms_individual / (ms_individual + (reps_per_ind - 1) * ms_error)
  
  list(pct_individual = pct_individual, pct_error = pct_error, f_stat = f_stat,
       df_individual = df_individual, df_error = df_error, pval = pval,
       repeatability = repeatability, aligned = aligned)
}

print_error_report <- function(label, result) {
  cat("\n", strrep("=", 60), "\n", label, "\n", strrep("=", 60), "\n", sep = "")
  cat(sprintf("Individual (biological) variance: %.2f%%\n", result$pct_individual))
  cat(sprintf("Measurement error variance:       %.2f%%\n", result$pct_error))
  cat(sprintf("F(%d,%d) = %.3f, permutation p = %.4f\n", result$df_individual, result$df_error, result$f_stat, result$pval))
  cat(sprintf("Repeatability index = %.3f\n", result$repeatability))
}

## ---- Per-landmark error breakdown -----------------------------------------
per_landmark_error <- function(coords_array, individual) {
  gpa_res <- gpagen(coords_array, print.progress = FALSE, curves = NULL)
  aligned <- gpa_res$coords
  k <- dim(aligned)[1]
  ind_levels <- unique(individual)
  landmark_error <- sapply(1:k, function(lm) {
    disp <- sapply(ind_levels, function(ind) {
      idx <- which(individual == ind)
      sum((aligned[lm, , idx[1]] - aligned[lm, , idx[2]])^2)
    })
    mean(disp)
  })
  data.frame(landmark = 1:k, mean_sq_replicate_disp = landmark_error)
}

## ============================================================================
## USAGE (requires ANTERIOR_REPLICATE.TPS / LATERAL_REPLICATE.TPS on disk):
##
## ant_input <- build_error_input("ANTERIOR__IMAGES_SAMPLE_AND_CONTROL.TPS", "ANTERIOR_REPLICATE.TPS")
## ant_result <- measurement_error_anova(ant_input$coords, ant_input$individual)
## print_error_report("ANTERIOR VIEW", ant_result)
## ant_lm_err <- per_landmark_error(ant_input$coords, ant_input$individual)
##
## lat_input <- build_error_input("LATERAL_IMAGES_4_SAMPLE_AND_CONTROL.TPS", "LATERAL_REPLICATE.TPS")
## lat_result <- measurement_error_anova(lat_input$coords, lat_input$individual)
## print_error_report("LATERAL VIEW", lat_result)
## lat_lm_err <- per_landmark_error(lat_input$coords, lat_input$individual)
## ============================================================================

if (file.exists("ANTERIOR_REPLICATE.TPS") && file.exists("LATERAL_REPLICATE.TPS")) {
  sink(tempfile())
    ant_input <- build_error_input("ANTERIOR__IMAGES_SAMPLE_AND_CONTROL.TPS", "ANTERIOR_REPLICATE.TPS")
    sink()
  ant_result <- measurement_error_anova(ant_input$coords, ant_input$individual)
  ant_lm_err <- per_landmark_error(ant_input$coords, ant_input$individual)
  write.csv(ant_lm_err[order(-ant_lm_err$mean_sq_replicate_disp), ],
            "anterior_per_landmark_error_geomorph.csv", row.names = FALSE)
  sink(tempfile())
  lat_input <- build_error_input("LATERAL_IMAGES_4_SAMPLE_AND_CONTROL.TPS", "LATERAL_REPLICATE.TPS")
  sink()
  lat_result <- measurement_error_anova(lat_input$coords, lat_input$individual)
  lat_lm_err <- per_landmark_error(lat_input$coords, lat_input$individual)
  write.csv(lat_lm_err[order(-lat_lm_err$mean_sq_replicate_disp), ],
            "lateral_per_landmark_error_geomorph.csv", row.names = FALSE)
  # Compile summary table dynamically using your exact function result outputs
  error_presentation_table <- tibble(
    `Anatomical Module Profile` = c("Anterior Symmetric Facial Plane", "Lateral Neurocranial Vault"),
    `Biological Variance (%)`   = c(round(ant_result$pct_individual, 1), round(lat_result$pct_individual, 1)),
    `Measurement Error Noise (%)` = c(round(ant_result$pct_error, 1), round(lat_result$pct_error, 1)),
    `F-value (df)`              = c(sprintf("%.2f (%d, %d)", ant_result$f_stat, ant_result$df_individual, ant_result$df_error), 
                                    sprintf("%.2f (%d, %d)", lat_result$f_stat, lat_result$df_individual, lat_result$df_error)),
    `P-value`                   = c(ifelse(ant_result$pval < 0.001, "< 0.001*", sprintf("%.4f*", ant_result$pval)),
                                    ifelse(lat_result$pval < 0.001, "< 0.001*", sprintf("%.4f*", lat_result$pval))),
    `Repeatability Index (R)`   = c(round(ant_result$repeatability, 3), round(lat_result$repeatability, 3))
  )
  
} else {
  cat("Replicate TPS files not found in the working directory.\n")
  cat("This script defines the functions only; see the USAGE block in the\n")
  cat("comments above once ANTERIOR_REPLICATE.TPS / LATERAL_REPLICATE.TPS are available.\n")
}
