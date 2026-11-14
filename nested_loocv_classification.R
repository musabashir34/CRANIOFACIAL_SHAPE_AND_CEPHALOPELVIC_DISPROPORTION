# nested_loocv_classification.R
# ---------------------------------------------------------------------------
# Leave-one-out cross-validated classification with NO information leakage from the left-out specimen.
#
# Within each leave-one-out fold:
#   1. PCA is refitted on the 159 training specimens only, and the left-out
#      specimen is projected onto those axes;
#   2. the classifier (linear discriminant / two-group CVA, equal priors) is
#      fitted on the training scores only;
#   3. (selected-k version) the number of PCs is chosen by an inner 10-fold
#      cross-validation inside the training set, with PCA refitted in each
#      inner fold as well.
# GPA is done once on the full sample: superimposition uses no group labels.
#
# Source AFTER the main pipeline (needs ant_sym_shape, lat_shape, ant_group,
# lat_group; geomorph loaded). Uses MASS (installed with R).
# Runtime: a few minutes.
# ---------------------------------------------------------------------------

set.seed(42)

nested_loocv <- function(shape_array, group, k_fixed = 10, k_grid = 2:20, inner_folds = 10) {
  X <- geomorph::two.d.array(shape_array)
  g <- factor(group, levels = c("C", "S"))
  n <- nrow(X)
  # restrict the search to PCs with real variance, so the reported k is the k used
  sd_full <- prcomp(X)$sdev
  k_grid  <- k_grid[k_grid <= sum(sd_full^2 / sum(sd_full^2) > 1e-8)]

  fit_predict <- function(Xtr, gtr, Xte, k) {
    pc  <- prcomp(Xtr, center = TRUE, scale. = FALSE)
    # keep only PCs with real variance: the anterior symmetric component has
    # 13 dimensions, so PCs 14+ are numerical noise that LDA cannot use
    n_real <- sum(pc$sdev^2 / sum(pc$sdev^2) > 1e-8)
    k   <- min(k, n_real)
    str <- pc$x[, 1:k, drop = FALSE]
    ste <- predict(pc, Xte)[, 1:k, drop = FALSE]
    fit <- MASS::lda(str, grouping = gtr, prior = c(0.5, 0.5))
    as.character(predict(fit, ste)$class)
  }

  pred_fixed <- pred_sel <- character(n)
  k_chosen   <- integer(n)

  for (i in seq_len(n)) {
    Xtr <- X[-i, , drop = FALSE]; gtr <- g[-i]

    # (a) fixed number of PCs, PCA refitted without specimen i
    pred_fixed[i] <- fit_predict(Xtr, gtr, X[i, , drop = FALSE], k_fixed)

    # (b) number of PCs selected by inner CV within the training set only
    folds <- sample(rep(seq_len(inner_folds), length.out = nrow(Xtr)))
    inner_acc <- sapply(k_grid, function(k) {
      hits <- unlist(lapply(seq_len(inner_folds), function(f) {
        te <- folds == f
        fit_predict(Xtr[!te, , drop = FALSE], gtr[!te], Xtr[te, , drop = FALSE], k) ==
          as.character(gtr[te])
      }))
      mean(hits)
    })
    k_chosen[i] <- k_grid[which.max(inner_acc)]
    pred_sel[i] <- fit_predict(Xtr, gtr, X[i, , drop = FALSE], k_chosen[i])
  }

  summarise <- function(pred) {
    truth <- as.character(g)
    c(accuracy    = 100 * mean(pred == truth),
      sensitivity = 100 * mean(pred[truth == "S"] == "S"),   # CPD = positive class
      specificity = 100 * mean(pred[truth == "C"] == "C"),
      correct     = sum(pred == truth))
  }

  list(fixed = summarise(pred_fixed), selected = summarise(pred_sel),
       k_median = median(k_chosen), k_range = range(k_chosen))
}

ant_ncv <- nested_loocv(ant_sym_shape, ant_group)
lat_ncv <- nested_loocv(lat_shape, lat_group)

nested_cv_results <- data.frame(
  view    = rep(c("Anterior", "Lateral"), each = 2),
  method  = rep(c("10 PCs, PCA refitted per fold", "PCs selected by inner CV per fold"), 2),
  n       = 160,
  correct = c(ant_ncv$fixed["correct"], ant_ncv$selected["correct"],
              lat_ncv$fixed["correct"], lat_ncv$selected["correct"]),
  accuracy    = round(c(ant_ncv$fixed["accuracy"], ant_ncv$selected["accuracy"],
                        lat_ncv$fixed["accuracy"], lat_ncv$selected["accuracy"]), 1),
  sensitivity = round(c(ant_ncv$fixed["sensitivity"], ant_ncv$selected["sensitivity"],
                        lat_ncv$fixed["sensitivity"], lat_ncv$selected["sensitivity"]), 1),
  specificity = round(c(ant_ncv$fixed["specificity"], ant_ncv$selected["specificity"],
                        lat_ncv$fixed["specificity"], lat_ncv$selected["specificity"]), 1),
  k_selected  = c(NA, sprintf("median %g (range %g-%g)", ant_ncv$k_median, ant_ncv$k_range[1], ant_ncv$k_range[2]),
                  NA, sprintf("median %g (range %g-%g)", lat_ncv$k_median, lat_ncv$k_range[1], lat_ncv$k_range[2]))
)

cat("\n=== Nested leave-one-out classification (no leakage) ===\n")
print(nested_cv_results)
cat("Original (PCA on full sample): anterior 68.8%, lateral 75.6%\n")
write.csv(nested_cv_results, "nested_loocv_classification.csv", row.names = FALSE)
