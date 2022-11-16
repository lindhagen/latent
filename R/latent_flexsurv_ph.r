#' @title Flexible spline proportional hazards model
#' @description Creates a model object for flexible spline proportional hazards
#' models.
#' @param formula An analysis formula. The dependent variable is assumed to be a
#' \link[survival]{Surv} object in data.
#' @param knots An array containing all knots, both boundary knots and inner
#' knots.
#' @param coef_names An optional character array of names for the regression
#' coefficients, excluding the spline coefficients. The latter are fixed as
#' \code{gamma0}, \code{gamma1} etc.
#' @return A \code{latent_flexsurv_ph} object.
#' @details The following generic functions can be invoked on the model object:
#' \itemize{
#'  \item{\link{latent_fit}}{   Fit the model.}
#'  \item{\link{latent_likelihood}}{   Compute likelihood.}
#'  \item{\link{latent_loglik_gradient_list}}{   Compute gradients of the
#'  log-likelihood.}
#'  \item{\link{latent_loglik_hessian_list}}{   Compute Hessians of the
#'  log-likelihood.}
#'  \item{\link{latent_simulate}}{   Simulate data from the model.}
#' }
#' @export
latent_flexsurv_ph <- function(formula, knots, coef_names = NULL) {
  nk <- length(knots) # Number of knots (internal + boundary).
  gamma_names <- sprintf("gamma%d", 0:(nk - 1))
  # No-intercept formula. Useful e.g. for linear predictors.
  if (length(all.vars(formula)) > 1) {
    formula0 <- formula %>% update(. ~ -1 + .)
  } else {
    # Formula has no covariates (e.g. y ~ 1). Replace it by (y ~ 1 - 1)
    # to avoid an intercept term in the model matrix.
    formula0 <- formula %>% update(. ~ -1)
  }
  # latent_flexsurv_ph object.
  obj <- list(formula = formula, formula0 = formula0, knots = knots,
    gamma_names = gamma_names, coef_names = coef_names)
  class(obj) <- "latent_flexsurv_ph"
  return (obj)
}

### ------------------------------ Generic functions ("public").
# Fit model.
latent_fit.latent_flexsurv_ph <- function(obj, data, weights) {
  ## Fit model.
  frm <- obj$formula
  knots <- obj$knots
  # Make the formula reside in current environment (for weights).
  environment(frm) <- environment()
  # Handle the issue  that flexssurvspline only accepts positive weights.
  # This also guarantees positive hazards at all data points, since the values
  # from the "wrong" subgroup are also present (although very weakly).
  weights <- pmax(weights, 1e-10)
  fit <- flexsurv::flexsurvspline(formula = frm,
    data = data,
    knots = knots %>% spline_inner,
    bknots = knots %>% spline_boundary,
    scale = "hazard",
    weights = weights)

  ## Extract results, and rename coefficients if requested.
  cof <- fit$coefficients
  if (!is.null(obj$coef_names)) {
    stopifnot(length(cof) == length(obj$gamma_names) + length(obj$coef_names))
    names(cof) <- c(obj$gamma_names, obj$coef_names)
  }

  # Return model parameters theta = c(gamma, beta) as a list.
  theta <- as.list(cof)
  return (theta)
}

# Likelihoods.
latent_likelihood.latent_flexsurv_ph <- function(obj, theta, data) {
  # Prepatations.
  prep <- latent_flexsurv_ph_prepare_diff(obj, theta, data)

  ### Likelihood.
  lik <- prep$dd %>%
    mutate(lik = exp(-cumhaz) * ifelse(y_evt, hazard, 1)) %>%
    pull(lik)

  return (lik)
}

# List of log-likelihood gradient vectors for the individuals in 'data'.
latent_loglik_gradient_list.latent_flexsurv_ph <- function(obj, theta, data) {
  # Preparations.
  prep <- latent_flexsurv_ph_prepare_diff(obj, theta, data)

  # Gradient list.
  grad_list <- lapply(1:nrow(data), function(i) {
    x <- prep$mm[i, ] %>% cbind
    b <- prep$spline_b[i, ] %>% cbind
    bp <- prep$spline_bp[i, ] %>% cbind
    gamma_bp <- sum(prep$gamma * bp) # Scalar product gamma * bp.
    H <- prep$dd$cumhaz[i]
    d <- prep$dd$y_evt[i]

    dloglik_dgamma <- (-H + d) * b + d / gamma_bp * bp
    dloglik_dbeta <- (-H + d) * x
    grad <- rbind(dloglik_dgamma, dloglik_dbeta) # Score vector.
    return (grad)
  })

  return (grad_list)
}

# List of Hessian matrices for the individuals in 'data'.
latent_loglik_hessian_list.latent_flexsurv_ph <- function(obj, theta, data) {
  # Preparations.
  prep <- latent_flexsurv_ph_prepare_diff(obj, theta, data)

  # Hessian list.
  hess_list <- lapply(1:nrow(data), function(i) {
    x <- prep$mm[i, ] %>% cbind
    b <- prep$spline_b[i, ] %>% cbind
    bp <- prep$spline_bp[i, ] %>% cbind
    xx <- tensor_square(x)
    bx <- tensor_prod(b, x)
    bb <- tensor_square(b)
    bpbp <- tensor_square(bp)
    gamma_bp <- sum(prep$gamma * bp) # Scalar product gamma * bp.
    H <- prep$dd$cumhaz[i]
    d <- prep$dd$y_evt[i]

    d2loglik_dgamma2 <- -H * bb - d / gamma_bp^2 * bpbp
    d2loglik_dgammadbeta <- -H * bx
    d2loglik_dbeta2 <- -H * xx
    # Build block matrix.
    hess_gamma <- cbind(d2loglik_dgamma2, d2loglik_dgammadbeta)
    hess_beta <- cbind(t(d2loglik_dgammadbeta), d2loglik_dbeta2)
    hess <- rbind(hess_gamma, hess_beta)
    return (hess)
  })
}

#' @title Simulate outcome
#' @description Simulates outcomes from a survival model,
#' see \link{latent_simulate} for details.
#' @export
latent_simulate.latent_flexsurv_ph <- function(obj, theta, data) {
  # Preparations.
  n <- nrow(data)
  knots <- obj$knots
  nk <- length(knots)
  theta %<>% unlist
  gamma <- theta[1:nk]
  beta <- theta[-(1:nk)]
  linpred = latent_linpred(obj$formula0, beta, data)

  # Matrix of gamma values.
  gamma_mat <- gamma %>%
    rep(times = n) %>%
    matrix(nrow = n, byrow = T)
  # PH: gamma0 is modified by adding the linear predictor.
  gamma_mat[, 1] %<>% add(linpred)

  # Simulate y.
  y_time_true <- flexsurv::rsurvspline(n = n,
    gamma = gamma_mat, knots = knots,
    scale = "hazard", timescale = "log")
  # Censor data randomly.
  # Since we don't know anything about the distribution of times, we generate
  # censoring times according to the same distribution as the outcome.
  # This gives 50% censoring.
  cns_time <- flexsurv::rsurvspline(n = n,
    gamma = gamma_mat, knots = knots,
    scale = "hazard", timescale = "log")
  y_time_obs <- pmin(y_time_true, cns_time)
  y_obs <- (y_time_true < cns_time)
  y <- Surv(y_time_obs, y_obs) # Censored data.

  return (y)

}

### --------------------------- Internal functions ("private").

# Prepare for differentiation (gradient, Hessian), by computing some
# auxiliary stuff.
latent_flexsurv_ph_prepare_diff <- function(obj, theta, data) {
  knots <- obj$knots
  nk <- length(knots)
  theta %<>% unlist
  gamma <- theta[1:nk]
  beta <- theta[-(1:nk)]
  yn <- all.vars(obj$formula)[1]

  # Tibble with individual-level stuff.
  dd <- tibble(
    # Outcome.
    y_time = data[[yn]][, "time"], # Follow-up time.
    y_evt = data[[yn]][, "status"], # Event indicator (TRUE/FALSE).
    log_y_time = log(y_time),

    # Hazard etc.
    linpred = latent_linpred(obj$formula0, beta, data), # Linear predictor.
    log_cumhaz = spline_fcn(log_y_time, gamma, knots) + linpred, # log(H).
    cumhaz = exp(log_cumhaz), # H.
    hazard = cumhaz * spline_fcn_deriv(log_y_time, gamma, knots)
  )
  # Safety check: the hazard should be positive.
  if (any(dd$hazard <= 0)) {
    stop("Negative hazard")
  }

  # x and spline-b matrices.
  mm <- model.matrix(obj$formula0, data = data)
  spline_b <- spline_basis_fcn_arr(dd$log_y_time, knots)
  spline_bp <- spline_basis_fcn_arr(dd$log_y_time, knots, deriv = T)

  # Return everything as a list.
  ret <- list(
    nk = nk, beta = beta, gamma = gamma,
    mm = mm, spline_b = spline_b, spline_bp = spline_bp,
    dd = dd
  )

  return (ret)
}
