% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/latent_glm.r
\name{latent_glm}
\alias{latent_glm}
\title{Generalized linear models}
\usage{
latent_glm(
  formula,
  family_name,
  link,
  coef_names = NULL,
  use_expected_info = F
)
}
\arguments{
\item{formula}{An analysis formula.}

\item{family_name}{A character variable defining the family.}

\item{link}{A character variable defining the link.}

\item{coef_names}{An optional character array of names for the regression
coefficients, including the intercept.}

\item{use_expected_info}{If \code{TRUE}, the expected Fischer information
will be used rather than the observed one when compuing Hessians.}
}
\value{
A \code{latent_glm} object.
}
\description{
Creates a model object for generalized linear models.
}
\details{
The following combinations of family and link are allowed:
\itemize{
 \item{\code{binomial}}{   Binomial outcome data. Links:
   \code{logit}, \code{probit}, \code{cauchit}, \code{cloglog}, \code{log},
   \code{identity}.}
 \item{\code{poisson}}{   Count outcome data. Links:
   \code{log}, \code{identity}, \code{sqrt}.}
}

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
