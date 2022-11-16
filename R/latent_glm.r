#' @title Generalized linear models
#' @description Creates a model object for generalized linear models.
#' @param formula An analysis formula.
#' @param family_name A character variable defining the family.
#' @param link A character variable defining the link.
#' @param coef_names An optional character array of names for the regression
#' coefficients, including the intercept.
#' @param use_expected_info If \code{TRUE}, the expected Fischer information
#' will be used rather than the observed one when compuing Hessians.
#' @return A \code{latent_glm} object.
#' @details The following combinations of family and link are allowed:
#' \itemize{
#'  \item{\code{binomial}}{   Binomial outcome data. Links:
#'    \code{logit}, \code{probit}, \code{cauchit}, \code{cloglog}, \code{log},
#'    \code{identity}.}
#'  \item{\code{poisson}}{   Count outcome data. Links:
#'    \code{log}, \code{identity}, \code{sqrt}.}
#' }
#'
#' The following generic functions can be invoked on the model object:
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
latent_glm <- function(formula, family_name, link,
  coef_names = NULL, use_expected_info = F) {

  if (family_name == "binomial") {

    family = binomial(link = link)
    ### New family functions for observed information.
    family$dvar_dmu <- function(mu) (1 - 2 * mu)
    ## Second derivetive of mu wrt eta.
    d2_fcn <- NULL # To be defined in below if statement.
    if (link == "logit") {
      d2_fcn <- function(x) -tanh(x/2) / (4 * cosh(x/2)^2)
    } else if (link == "probit") {
      d2_fcn <- function(x) -x * dnorm(x)
    } else if (link == "cauchit") {
      d2_fcn <- function(x) -2 * x / (pi * (1 + x^2)^2)
    } else if (link == "cloglog") {
      d2_fcn <- function(x) exp(-exp(x)) * exp(x) * (1 - exp(x))
    } else if (link == "log") {
      d2_fcn <- function(x) exp(x)
    } else if (link == "identity") {
      d2_fcn <- function(x) rep(0, length(x))
    } else {
      stop("Unexpected link") # Already have an error from binomial().
    }
  } else if (family_name == "poisson") {

    family = poisson(link = link)
    ### New family functions for observed information.
    family$dvar_dmu <- function(mu) rep(1, length(mu))
    ## Second derivetive of mu wrt eta.
    d2_fcn <- NULL # To be defined in below if statement.
    if (link == "log") {
      d2_fcn <- function(x) exp(x)
    } else if (link == "identity") {
      d2_fcn <- function(x) rep(0, length(x))
    } else if (link == "sqrt") {
      d2_fcn <- function(x) rep(2, length(x))
    } else {
      stop("Unexpected link") # Already have an error from poisson().
    }
  } else {
    stop("Illegal latent_glm family")
  }
  family$d2mu_deta2 <- d2_fcn

  # latent_glm object.
  obj <- list(
    formula = formula,
    family_name = family_name,
    family = family,
    coef_names = coef_names,
    use_expected_info = use_expected_info)
  class(obj) <- "latent_glm"
  return (obj)
}

### ------------------------------ Generic functions ("public").
# Fit model.
latent_fit.latent_glm <- function(obj, data, weights) {
  ## Fit model.
  frm <- obj$formula
  # Make the formula reside in current environment (for weights).
  environment(frm) <- environment()
  # Suppress warning about non-integer weights.
  fit <- suppressWarnings(glm(
    formula = frm, data = data, family = obj$family, weights = weights))

  ## Extract results, and rename coefficients if requested.
  cof <- fit$coef
  if (!is.null(obj$coef_names)) {
    stopifnot(length(obj$coef_names) == length(cof))
    names(cof) <- obj$coef_names
  }

  # Return model parameters theta (= beta; there are no additional
  # model specifics) as a list.
  theta <- as.list(cof)
  return (theta)
}

# Likelihoods.
latent_likelihood.latent_glm <- function(obj, theta, data) {
  # Don't know the name of the outcome variable. Extract it from the formula.
  yn <- all.vars(obj$formula)[1]
  y <- data[[yn]]
  # Linear predictor.
  eta <- latent_linpred(formula = obj$formula, beta = theta, data = data)
  mu <- obj$family$linkinv(eta) # g(mu) = eta.

  if (obj$family_name == "binomial") {
    lik <- ifelse(y, mu, 1 - mu)
  } else if (obj$family_name == "poisson") {
    lik <- dpois(y, lambda = mu)
  } else {
    stop("Unknown latent_glm family")
  }

  return (lik)
}

# List of log-likelihood gradient vectors for the individuals in 'data'.
latent_loglik_gradient_list.latent_glm <- function(obj, theta, data) {
  mm <- model.matrix(obj$formula, data = data)
  dd <- latent_glm_prepare_diff(obj, theta, data)

  grad_list <- lapply(1:nrow(data), function(i) {
    x <- mm[i, ] %>% cbind
    y <- dd$y[i]
    mu <- dd$mu[i]
    v <- dd$v[i]
    dmu_deta = dd$dmu_deta[i]
    dloglik_dbeta <- (y - mu) / v * dmu_deta * x
    grad <- rbind(dloglik_dbeta) # Score vector.
    return (grad)
  })
  return (grad_list)
}

# List of Hessian matrices for the individuals in 'data'.
latent_loglik_hessian_list.latent_glm <- function(obj, theta, data) {
  mm <- model.matrix(obj$formula, data = data)
  dd <- latent_glm_prepare_diff(obj, theta, data)

  hess_list <- lapply(1:nrow(data), function(i) {
    x <- mm[i, ] %>% cbind
    y <- dd$y[i]
    mu <- dd$mu[i]
    v <- dd$v[i]
    vp <- dd$dv_dmu[i]
    dmu_deta = dd$dmu_deta[i]
    d2mu_deta2 = dd$d2mu_deta2[i]

    # Scalar part of d2(loglik)/d(theta)2, i.e. without x * x.
    scalar <- (-1 / v) * (dmu_deta)^2
    if (!obj$use_expected_info) {
      scalar %<>% add(
        ((y - mu) / v) * (d2mu_deta2 - (vp / v) * (dmu_deta)^2))
    }
    d2loglik_dbeta2 <- scalar * tensor_square(x)

    hess <- d2loglik_dbeta2
    return (hess)
  })

  return (hess_list)
}


#' @title Simulate outcome
#' @description Simulates outcomes from a generalized linear model,
#' see \link{latent_simulate} for details.
#' @export
latent_simulate.latent_glm <- function(obj, theta, data) {
  eta <- latent_linpred(formula = obj$formula, beta = theta, data = data)
  mu <- obj$family$linkinv(eta)

  if (obj$family_name == "binomial") {
    y <- rbern(mu)
  } else if (obj$family_name == "poisson") {
    y <- rpois(n = length(mu), lambda = mu)
  } else {
    stop("Unknown latent_glm family")
  }

  return (y)
}

### --------------------------- Internal functions ("private").

# Prepare for differentiation (gradient, Hessian), by computing some
# auxiliary stuff.
latent_glm_prepare_diff <- function(obj, theta, data) {
  yn <- all.vars(obj$formula)[1]
  fam <- obj$family
  ret <- data %>% mutate(
    y = data[[yn]],
    # Linear predictor.
    eta = latent_linpred(formula = obj$formula, beta = theta, data = data),
    mu = fam$linkinv(eta),
    v = fam$variance(mu), # Variance = v(mu).
    dv_dmu = fam$dvar_dmu(mu), # v'(mu).
    dmu_deta = fam$mu.eta(eta),
    d2mu_deta2 = fam$d2mu_deta2(eta)
  )
  return (ret)
}

