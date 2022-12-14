% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/latent_linear.r
\name{latent_linear}
\alias{latent_linear}
\title{Linear model}
\usage{
latent_linear(formula, coef_names = NULL)
}
\arguments{
\item{formula}{An analysis formula.}

\item{coef_names}{An optional character array of names for the regression
coefficients, including the intercept.}
}
\value{
A \code{latent_linear} object.
}
\description{
Created a model object for linear models.
}
\details{
The following generic functions can be invoked on the model object:
\itemize{
 \item{\link{latent_fit}}{   Fit the model.}
 \item{\link{latent_likelihood}}{   Compute likelihood.}
 \item{\link{latent_loglik_gradient_list}}{   Compute gradients of the
 log-likelihood.}
 \item{\link{latent_loglik_hessian_list}}{   Compute Hessians of the
 log-likelihood.}
 \item{\link{latent_simulate}}{   Simulate data from the model.}
}
}
