#' @title Ordinal proportional odds models
#' @description Creates a model object for ordinal models, see the package
#' vignette for details. The outcome variable should be a factor with 3+ levels.
#' It need not be \link{ordered} in R, the order of the levels defined the
#' ordinality.
#' @param formula An analysis formula.
#' @param K The number of levels (at least 3) in the ordinal response variable.
#' @param coef_names An optional character array of names for the regression
#' coefficients, excluding the intercepts. The latter are fixed as
#' \code{alpha1}, \code{alpha2} etc.
#' @return A \code{latent_ordinal} object.
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
latent_ordinal <- function(formula, K, coef_names = NULL) {
  stopifnot (K >= 3)
  Km1 <- K - 1
  alpha_names <- sprintf("alpha%d", 1:Km1)
  # No-intercept formula. Useful e.g. for linear predictors.
  if (length(all.vars(formula)) > 1) {
    formula0 <- formula %>% update(. ~ -1 + .)
  } else {
    # Formula has no covariates (e.g. y ~ 1). Replace it by (y ~ 1 - 1)
    # to avoid an intercept term in the model matrix.
    formula0 <- formula %>% update(. ~ -1)
  }
  # latent_ordinal object.
  obj <- list(formula = formula, formula0 = formula0, K = K,
    alpha_names = alpha_names, coef_names = coef_names)
  class(obj) <- "latent_ordinal"
  return (obj)
}

### ------------------------------ Generic functions ("public").
# Fit model.
latent_fit.latent_ordinal <- function(obj, data, weights) {
  ## Fit model.
  frm <- obj$formula
  # Make the formula reside in current environment (for weights).
  environment(frm) <- environment()

  fit <- suppressWarnings( # polr warns about non-integer weights.
    MASS::polr(formula = frm,
      data = data,
      weights = weights)
  )

  ## Extract results, and rename coefficients if requested.
  # polr has a different sign convention, need to change sign of alpha.
  alpha <- -fit$zeta
  beta <- fit$coefficients
  names(alpha) <- obj$alpha_names
  if (!is.null(obj$coef_names)) {
    stopifnot(length(beta) == length(obj$coef_names))
    names(beta) <- obj$coef_names
  }
  cof <- c(alpha, beta)

  # Return model parameters theta = c(gamma, beta) as a list.
  theta <- as.list(cof)
  return (theta)
}

# Likelihoods.
latent_likelihood.latent_ordinal <- function(obj, theta, data) {
  # Prepatations.
  prep <- latent_ordinal_prepare_diff(obj, theta, data)

  # Extract likelihood (simple!).
  dd <- prep$dd
  lik <- dd$lik

  return (lik)
}

# List of log-likelihood gradient vectors for the individuals in 'data'.
latent_loglik_gradient_list.latent_ordinal <- function(obj, theta, data) {
  # Preparations.
  prep <- latent_ordinal_prepare_diff(obj, theta, data)
  mm <- prep$mm
  dd <- prep$dd
  en <- prep$en
  An <- prep$An
  K <- prep$K
  lik <- dd$lik # Likelihood.

  # Gradient list.
  grad_list <- lapply(1:nrow(data), function(i) {
    x <- mm[i, ] %>% cbind
    # epsilon and A column vectors (length K-1).
    eps <- dd[i, en[-K]] %>% as.matrix %>% t
    A <- dd[i, An] %>% as.matrix %>% t

    # Derivatives of the likelihood.
    dlik_dalpha <- eps * A
    dlik_dbeta <- sum(eps * A) * x
    dlik_dtheta <- rbind(dlik_dalpha, dlik_dbeta)
    # Convert to derivative of the log-likelihood (score vector).
    dloglik_dtheta <- logdiff_grad_plain2log(Fx = lik[i], grad_F = dlik_dtheta)

    return (dloglik_dtheta)
  })

  return (grad_list)
}

# List of Hessian matrices for the individuals in 'data'.
latent_loglik_hessian_list.latent_ordinal <- function(obj, theta, data) {
  # Preparations.
  prep <- latent_ordinal_prepare_diff(obj, theta, data)
  mm <- prep$mm
  dd <- prep$dd
  en <- prep$en
  Bn <- prep$Bn
  K <- prep$K
  # Need likelihood and gradients for logarithmic differentiation.
  lik <- dd$lik # Likelihood.
  dloglik_dtheta_list <- latent_loglik_gradient_list(obj, theta, data)

  # Hessian list.
  hess_list <- lapply(1:nrow(data), function(i) {
    x <- mm[i, ] %>% cbind
    # epsilon and C as plain vectors (length K-1).
    eps <- dd[i, en[-K]] %>% as.numeric
    B <- dd[i, Bn] %>% as.numeric

    epsB <- eps * B # Element-wise products, length K-1.
    # Derivatives of the likelihood.
    d2lik_dalpha2 <- diag(epsB)
    d2lik_dbetadalpha <- tensor_prod(x, epsB)
    d2lik_dbeta2 <- sum(epsB) * tensor_square(x)
    # Build block matrix.
    upper <- cbind(d2lik_dalpha2, t(d2lik_dbetadalpha))
    lower <- cbind(d2lik_dbetadalpha, d2lik_dbeta2)
    d2lik_dtheta2 <- rbind(upper, lower)
    # Convert to derivative of the log-likelihood.
    d2loglik_dtheta2 <- logdiff_hess_plain2log_hybrid(
      Fx = lik[i],
      grad_logF = dloglik_dtheta_list[[i]],
      hess_F = d2lik_dtheta2)

    return (d2loglik_dtheta2)
  })
}

#' @title Simulate outcome
#' @description Simulates outcomes from an ordinal model,
#' see \link{latent_simulate} for details.
#' @export
latent_simulate.latent_ordinal <- function(obj, theta, data) {
  # Preparations.
  n <- nrow(data)
  K <- obj$K
  Km1 <- K - 1
  theta %<>% unlist
  alpha <- theta[1:Km1]
  beta <- theta[-(1:Km1)]
  linpred = latent_linpred(obj$formula0, beta, data)

  u <- runif(n)
  y <- rep(NA, n)
  for (k in 1:Km1) {
    # Model probabilities, P(y <= k).
    p.le.k = 1 / (1 + exp(alpha[k] + linpred))
    # Outcome found if this is the first time that u < p.
    y <- ifelse(is.na(y) & (u < p.le.k), k, y)
  }
  # If y is still NA, no outcome ws found among 1, 2, ..., K-1, so the
  # outcome has to be K.
  y <- ifelse(is.na(y), K, y)

  # Convert to ordered factor.
  y %<>% factor(ordered = T)

  return (y)
}

### --------------------------- Internal functions ("private").

# Prepare for differentiation (gradient, Hessian), by computing some
# auxiliary stuff.
latent_ordinal_prepare_diff <- function(obj, theta, data) {
  K <- obj$K
  Km1 <- K - 1
  theta %<>% unlist
  alpha <- theta[1:Km1]
  beta <- theta[-(1:Km1)]
  yn <- all.vars(obj$formula)[1]

  ## Tibble with individual-level stuff.
  y <- as.numeric(data[[yn]]) # The outcome as a numerical variable (1 - K).
  dd <- tibble(y = y)
  en <- sprintf("eps%d", 1:K)
  An <- sprintf("A%d", 1:Km1)
  Bn <- sprintf("B%d", 1:Km1)
  Cn <- sprintf("C%d", 1:K)
  linpred <- latent_linpred(formula = obj$formula0, beta, data)
  for (k in 1:K) {
    eps <- case_when(
      k == y ~ 1,
      k == y - 1 ~ -1,
      TRUE ~ 0
    )
    dd[[en[k]]] <- eps
    if (k < K) {
      nu <- (alpha[k] + linpred) / 2
      C <- 1 / (1 + exp(2 * nu))
      A <- -1 / (4 * cosh(nu)^2)
      B <- -A * tanh(nu)
      # Write back to tibble.
      dd[[An[k]]] <- A
      dd[[Bn[k]]] <- B
      dd[[Cn[k]]] <- C
    }
  }
  dd[[Cn[K]]] <- 1

  # Add likelihood.
  dd %<>% mutate(lik = 0)
  for (k in 1:K) {
    enk <- en[k]
    Cnk <- Cn[k]
    dd$lik %<>% add(dd[[enk]] * dd[[Cnk]])
  }

  # Sort columns.
  dd %<>% select(y, starts_with("eps"), starts_with("C"), starts_with("A"),
    starts_with("B"), lik)

  # Model matrix.
  mm <- model.matrix(obj$formula0, data = data)

  # Return everything as a list.
  ret <- list(
    K = K, alpha = alpha, beta = beta,
    en = en, An = An, Bn = Bn, Cn = Cn,
    mm = mm,
    dd = dd
  )

  return (ret)
}
