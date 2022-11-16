# Extract linear predictors given analysis formula, data, and regression
# coefficients. By "regression coefficients", we mean the linear component of
# the model (beta), not other model parameters (like sigma in a linear model).
# However, they are still supposed to come as list.
latent_linpred <- function(formula, beta, data) {
  if (length(beta) > 0) {
    # Remove dependent variable (y) from the formula, since we want to be able
    # to do without it.
    frm <- update(formula, NULL ~ .)
    mm <- model.matrix(frm, data = data)
    beta_c <- beta %>% unlist %>% cbind
    stopifnot(ncol(mm) == nrow(beta_c)) # Dimension consistency check.

    lin_pred <- mm %*% beta_c # Matrix multiplication.
    lin_pred <- lin_pred[, 1] # Un-matrix.
  } else {
    # No regression coefficients, don't add anything anywhere.
    lin_pred <- rep(0, nrow(data))
  }
  return (lin_pred)
}
