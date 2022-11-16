#' @title Fit a model
#' @description Fits a model to data using weighted maximum likelihood.
#' @param obj A model object.
#' @param data A data frame.
#' @param weights An array of non-negative weights of the same length as the
#' number of rows in \code{data}.
#' @return A list with fitted parameters.
#' @export
latent_fit <- function(obj, ...) {
  UseMethod("latent_fit", obj)
}

#' @title Compute likelihood
#' @description Computes the likelihood for a model object on given data.
#' @param obj A model object.
#' @param theta A list of model parameters.
#' @param data A data frame.
#' @return An array with individual likelihoods.
#' @export
latent_likelihood <- function(obj, ...) {
  UseMethod("latent_likelihood", obj)
}

#' @title Compute likelihood gradients
#' @description Computes the gradient (vector with first-order derivatives)
#' of the log-likelihood for a model object on given data.
#' @param obj A model object.
#' @param theta A list of model parameters.
#' @param data A data frame.
#' @return A list with  gradients for the given individuals.
#' @export
latent_loglik_gradient_list <- function(obj, ...) {
  UseMethod("latent_loglik_gradient_list", obj)
}

#' @title Compute likelihood Hessians
#' @description Computes the Hessian (matrix with second-order derivatives)
#' of the log-likelihood for a model object on given data.
#' @param obj A model object.
#' @param theta A list of model parameters.
#' @param data A data frame.
#' @return A list with  Hessians for the given individuals.
#' @export
latent_loglik_hessian_list <- function(obj, ...) {
  UseMethod("latent_loglik_hessian_list", obj)
}

# Simulate outcome.
# Expected parameters:
#  - theta = Parameter vector.
#  - data = Data set with covariates, but without outcome.
# Should return an array of simulated outcomes.
#' @title Simulate outcome from a model
#' @description Simulates the outcome \code{y} for a model with given
#' @param obj A model object.
#' @param theta A list of model parameters.
#' @param data A data frame with covariates for the simulation.
#' @return An array with simulated outcomes.
#' @export
latent_simulate <- function(obj, ...) {
  UseMethod("latent_simulate", obj)
}
