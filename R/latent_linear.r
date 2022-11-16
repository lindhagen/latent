### Linear model.

#' @title Linear model
#' @description Created a model object for linear models.
#' @param formula An analysis formula.
#' @param coef_names An optional character array of names for the regression
#' coefficients, including the intercept.
#' @return A \code{latent_linear} object.
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
latent_linear <- function(formula, coef_names = NULL) {
  obj <- list(formula = formula, coef_names = coef_names)
  class(obj) <- "latent_linear"
  return (obj)
}

### ------------------------------ Generic functions ("public").

# Fit linear model.
latent_fit.latent_linear <- function(obj, data, weights) {
  ## Fit model.
  frm <- obj$formula
  # Make the formula reside in current environment. O/w, lm won't recognize
  # 'weights'. This may not be state of the art, but it seems to work!
  environment(frm) <- environment()
  fit <- lm(formula = frm, data = data, weights = weights)

  ## Extract results, and rename coefficients if requested.
  cof <- fit$coef
  if (!is.null(obj$coef_names)) {
    stopifnot(length(obj$coef_names) == length(cof))
    names(cof) <- obj$coef_names
  }
  # Residual variance. Due to weighting, we can't extract this directly from
  # the lm object. Do it "by hand" instead.
  resid <- residuals(fit)
  sig2 <- weighted.mean(resid^2, w = weights)

  # Return model parameters theta = (beta, sigma2) as a list.
  theta <- c(
    as.list(cof),
    sigma2 = sig2 # Residual variance always called 'sigma2'.
  )
  return (theta)
}

# Likelihoods.
latent_likelihood.latent_linear <- function(obj, theta, data) {
  pred <- latent_linear_predict(obj, theta, data)
  lik <- dnorm(data$y, mean = pred, sd = sqrt(theta$sigma2))
  return (lik)
}

# List of log-likelihood gradient vectors for the individuals in 'data'.
latent_loglik_gradient_list.latent_linear <- function(obj, theta, data) {
  s2 <- theta$sigma2 # Short-hand.
  mm <- model.matrix(obj$formula, data = data)
  dd <- data %>%
    mutate(resid = latent_linear_resid(obj, theta, data))
  grad.list <- lapply(1:nrow(data), function(i) {
    x <- mm[i, ] %>% cbind
    r <- dd$resid[i]
    dloglik.dbeta <- (r / s2) * x
    dloglik.dsig2 <- (r^2 - s2) / (2 * s2^2)
    grad <- rbind(dloglik.dbeta, sigma2 = dloglik.dsig2) # Gradient vector.
    return (grad)
  })
  return (grad.list)
}

# List of Hessian matrices for the individuals in 'data'.
latent_loglik_hessian_list.latent_linear <- function(obj, theta, data) {
  s2 <- theta$sigma2 # Short-hand.
  mm <- model.matrix(obj$formula, data = data)
  dd <- data %>%
    mutate(resid = latent_linear_resid(obj, theta, data))
  hess.list <- lapply(1:nrow(data), function(i) {
    x <- mm[i, ] %>% cbind
    xx <- tensor_square(x)
    r <- dd$resid[i]
    d2loglik.dbeta2 <- -xx / s2
    d2loglik.dsig2.dbeta <- (-r / s2^2) * x
    d2loglik.dsig22 <- (s2 - 2* r^2) / (2 * s2^3)
    colnames(d2loglik.dsig2.dbeta) <- "sigma2"
    hess.up <- cbind(d2loglik.dbeta2, d2loglik.dsig2.dbeta)
    hess.lo <- cbind(t(d2loglik.dsig2.dbeta), d2loglik.dsig22)
    hess <- rbind(hess.up, hess.lo)
    return (hess)
  })
  return (hess.list)
}

#' @title Simulate outcome
#' @description Simulates outcomes from a linear model,
#' see \link{latent_simulate} for details.
#' @export
latent_simulate.latent_linear <- function(obj, theta, data) {
  y.sim <- latent_linear_predict(obj, theta, data) +
    rnorm(nrow(data), sd = sqrt(theta$sigma2))
  return (y.sim)
}



### --------------------------- Internal functions ("private").

# Predictions from linear model.
# Simple wrapper around general linear predictor function, removing sigma.
latent_linear_predict <- function(obj, theta, data) {
  beta <- theta
  beta$sigma2 <- NULL #  Only regression coefficients.
  pred <- latent_linpred(formula = obj$formula, beta = beta, data = data)
  return (pred)
}

# Residuals from linear model.
latent_linear_resid <- function(obj, theta, data) {
  pred <- latent_linear_predict(obj, theta, data)
  resid <- data$y - pred
  return (resid)
}
