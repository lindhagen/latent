#' @title Effect estimation in latent subgroups
#' @description Main function for performing effect estimation in latent
#' subgroups.
#' @param data A data frame containing the data.
#' @param modelS A model object  (e.g. \link{latent_glm}) for subgroup
#' # membership (dichotomous outcome data).
#' @param model0 Outcome model for individuals not belonging to the latent
#' subgroup (untreatable patients).
#' @param model1 Outcome model for individuals belonging to the latent
#' subgroup (treatable patients).
#' @param do_em If \code{TRUE}, the EM algorithm will be used. Otherwise,
#' the likelihood is maximized using the general-purpose optimization function
#' \link{optim} with \code{method = BFGS}.
#' @param max_em_iter Maximum number of iterations in the EM algorithm. If the
#' algorithm doesn't converge within this number of iterations, an error is
#' raised.
#' @param loglik_conv_tol Convergence tolerance criterion for the EM algorithm.
#' If the log-likelohood increases less than this amount between two successive
#' iterations, the algorithm is considered to have converged.
#' @return A list with two fields:
#' \itemize{
#'  \item{\code{theta}}{   An array \eqn{\theta = (\theta_S, \theta_0, \theta_1)}
#'  with fitted parameters from all three models.}
#'  \item{\code{info}}{   The Fisher information matrix, i.e. the negavive
#'  Hessian of the log-likelihood with respect to the full set of model
#'  parameters \eqn{\theta} (all three models).}
#' }
#' @details Dichotomous variables should be 0/1 coded. For the latent subgroup,
#' 0 will correspond to untreatable patients, and 1 to treatable ones.
#'
#' To avoid name conflicts, the variable names in \code{data} should
#' not start with a dot (\code{.}). The reason is that this function will add
#' some variables for internal use, all starting with a dot.
#'
#' For further details on the usage of this function, see the package vignette.
#' @importFrom magrittr %>% %<>% add
#' @importFrom dplyr tibble mutate select pull case_when
#' @export
latent_main <- function(data,
  modelS,
  model0,
  model1,
  do_em = T,
  max_em_iter = 250,
  loglik_conv_tol = 1e-5 # log-likelihood convergence criterion for EM.
) {

  ### Utility function: Enrich the given data set with pi/likelihood stuff.
  ### Model objects are global are variables.
  setup.likelihood.data <- function(data, theta) {
    ret <- data %>%
      mutate(
        # Markers for possibility of subgroup membership.
        # delta = 1 unless g is known and opposite.
        delta0 = ifelse(!is.na(.grp) & (.grp == 1), 0, 1),
        delta1 = ifelse(!is.na(.grp) & (.grp == 0), 0, 1),
        ## Model likelihoods.
        # Likelihood of g = 0, 1.
        MS0 = latent_likelihood(modelS, theta = theta$thetaS,
          data = data_g0),
        MS1 = 1 - MS0,
        # Outcome models.
        M0 = latent_likelihood(model0, theta = theta$theta0, data = data),
        M1 = latent_likelihood(model1, theta = theta$theta1, data = data),
        lik = delta0 * MS0 * M0 + delta1 * MS1 * M1,
        loglik = log(lik)
      )
    return (ret)
  }

  ### We shall repeatedly need to access the subgroup variable. For covenience,
  ### we create a copy of it, called '.grp'.
  xn_subgroup <- all.vars(modelS$formula)[1] # Subgroup variable name.
  data$.grp = data[[xn_subgroup]] # The copy.
  # Also, create copies of the data set, where the subgroup is identically
  # zero or one. This will be useful for certain likelihoods.
  data_g0 <- data
  data_g0[[xn_subgroup]] <- 0
  data_g1 <- data
  data_g1[[xn_subgroup]] <- 1

  ### Find start values (EM or numerical optimization).
  # For thetaS, we fit the model on treatable patients, where the latent
  # subgroup is observed.
  #
  # For theta0 and theta1, we fit weighted models, where the weights for
  # unobserved subgroup are based on the predicted probabilities from the
  # subgroup model described above.
  thetaS_start <- latent_fit(modelS, data = data, weights = rep(1, nrow(data)))
  # Compute predicted subgroup probabilities.
  dd_start <- data_g1
  dd_start %<>% mutate(
    .pg = latent_likelihood(modelS, theta = thetaS_start, data = dd_start),
    # Weight = subgroup (0/1) when observed.
    .wt = ifelse(is.na(.grp), .pg, .grp),
    # Keep the weights away from 0 and 1 (e.g. to prevent negative hazards in
    # survival spline models).
    .wt = pmax(.wt, 0.001),
    .wt = pmin(.wt, 0.999)
  )
  theta0_start <- latent_fit(model0, data = dd_start, weights = 1 - dd_start$.wt)
  theta1_start <- latent_fit(model1, data = dd_start, weights = dd_start$.wt)
  theta_start <- list(
    thetaS = thetaS_start,
    theta0 = theta0_start,
    theta1 = theta1_start)
  theta_start_flat <- flatten.listlist(theta_start)

  # Convenience text array 'si.names' = nice names for gradients/Hessians
  dS <- length(theta_start$thetaS)
  d0 <- length(theta_start$theta0)
  d1 <- length(theta_start$theta1)
  si_names <- rep("", dS + d0 + d1)
  si_names[1] <- "thetaS"
  si_names[1 + dS] <- "theta0"
  si_names[1 + dS + d0] <- "theta1"

  # Either EM + closed form standard errors, or numerical methods.
  if (do_em) {
    ### EM algorithm.
    # EM start-up.
    theta_old <- theta_start
    theta_old_flat <- theta_start_flat
    loglik_old <- setup.likelihood.data(data = data, theta = theta_old) %>%
      pull(loglik) %>%
      sum

    # Data to track EM iterations.
    dd.em <- tibble(dummy = rep(NA, max_em_iter))
    for (xn in names(theta_start_flat)) {
      dd.em[[xn]] <- NA
      dd.em[[xn]][1] <- theta_start_flat[[xn]]
    }
    dd.em %<>% select(-dummy)
    em_convergence <- FALSE # No convergence yet.

    # Main EM loop. It continues until convergence or max number of iterations
    # (failure).
    for (em_ix in 2:max_em_iter) {
      ## E step.
      # The job here is to compute the posterior probabilities 'w' of belonging to
      # group 1 using Bayes's theorem. To be used as weights in regression models
      # in the M step.
      #
      # Model likelihoods.
      MS0 <- latent_likelihood(modelS, theta = theta_old$thetaS,
        data = data_g0)
      MS1 <- 1 - MS0
      M0 <- latent_likelihood(model0, theta = theta_old$theta0, data = data)
      M1 <- latent_likelihood(model1, theta = theta_old$theta1, data = data)
      # Write lots of probability info to the working data set 'dd'.
      dd <- data %>%
        mutate(
          # Prior probabilities.
          pi0 = ifelse(is.na(.grp), MS0, 1 - .grp),
          pi1 = ifelse(is.na(.grp), MS1, .grp),
          # Posterios probabilites.
          w0t = pi0 * M0, # w-tilde (unnormalized posterior prob's).
          w1t = pi1 * M1,
          # Posterior probability of belonging to group 1.
          w = w1t / (w0t + w1t))

      ## M step.
      ## Three weighted regression models for thetaS, theta0 and theta1,
      ## using the weigths from the E step.
      # Subgroup model.
      # Augment the data so that each patient contributes two rows, with
      # (subgroup, weight) = (0, 1 - w) and (1, w), respectively.
      ddS0 <- dd %>% mutate(wt = 1 - w)
      ddS0[[xn_subgroup]] <- 0
      ddS1 <- dd %>% mutate(wt = w)
      ddS1[[xn_subgroup]] <- 1
      ddS <- rbind(ddS0, ddS1)
      thetaS <- latent_fit(modelS, data = ddS, weights = ddS$wt)
      # Outcome models.
      theta0 <- latent_fit(model0, data = data, weights = 1 - dd$w)
      theta1 <- latent_fit(model1, data = data, weights = dd$w)
      # Summarize in theta list.
      theta_new <- list(
        thetaS = thetaS,
        theta0 = theta0,
        theta1 = theta1
      )
      theta_new_flat <- flatten.listlist(theta_new)
      stopifnot(identical(names(theta_new_flat), names(theta_old_flat)))

      # Write back to result tibble and check convergence.
      dd.em[em_ix, ] <- as.list(theta_new_flat)
      loglik_new <- setup.likelihood.data(data = data, theta = theta_new) %>%
        pull(loglik) %>%
        sum
      if (abs(loglik_new - loglik_old) < loglik_conv_tol) {
        em_convergence <- TRUE
        break # Exit EM loop.
      }
      # No convergence yet, keep looping.
      theta_old <- theta_new
      theta_old_flat <- theta_new_flat
      loglik_old <- loglik_new
    } # End EM loop.
    if (!em_convergence) {
      # Convergence failed. Throw error and abort.
      stop(sprintf("EM algorithm didn't converge within %d iterations\n.",
        max_em_iter))
    }
    # Output from EM algorithm.
    theta_em <- theta_new
    theta_em_flat <- theta_new_flat


    ## Closed-form first (gradient) and second (Hessian) derivatives of the
    ## log-likelihood. We call these _em, alghough they have nothing to do
    ## with the EM algorithm, they just go together.
    ## Individual results are stored as lists, with one item per individual.
    ## Eventually, we will sum over these (Reduce).

    ## Individual derivatives of the model-specific log-likelihood.
    ## Stored as lists with one vector/matrix per individual.

    # Likelihoods.
    dd_lik <- setup.likelihood.data(data = data, theta = theta_em)

    ## Log-likelihood derivatives.
    # Gradients.
    loglik_gradS0_list <- latent_loglik_gradient_list(
      modelS, theta = theta_em$thetaS, data = data_g0)
    loglik_gradS1_list <- latent_loglik_gradient_list(
      modelS, theta = theta_em$thetaS, data = data_g1)
    loglik_grad0_list <- latent_loglik_gradient_list(
      model0, theta = theta_em$theta0, data = data)
    loglik_grad1_list <- latent_loglik_gradient_list(
      model1, theta = theta_em$theta1, data = data)
    # Hessians.
    loglik_hessS0_list <- latent_loglik_hessian_list(
      modelS, theta = theta_em$thetaS, data = data_g0)
    loglik_hessS1_list <- latent_loglik_hessian_list(
      modelS, theta = theta_em$thetaS, data = data_g1)
    loglik_hess0_list <- latent_loglik_hessian_list(
      model0, theta = theta_em$theta0, data = data)
    loglik_hess1_list <- latent_loglik_hessian_list(
      model1, theta = theta_em$theta1, data = data)

    ## Likelihood derivatives.
    # Gradients.
    lik_gradS0.list <- map2(as.list(dd_lik$MS0), loglik_gradS0_list,
      logdiff_grad_log2plain)
    lik_gradS1.list <- map2(as.list(dd_lik$MS1), loglik_gradS1_list,
      logdiff_grad_log2plain)
    lik_grad0.list <- map2(as.list(dd_lik$M0), loglik_grad0_list,
      logdiff_grad_log2plain)
    lik_grad1.list <- map2(as.list(dd_lik$M1), loglik_grad1_list,
      logdiff_grad_log2plain)
    # Hessians.
    lik_hessS0.list <- map3(as.list(dd_lik$MS0),
      loglik_gradS0_list, loglik_hessS0_list,
      logdiff_hess_log2plain)
    lik_hessS1.list <- map3(as.list(dd_lik$MS1),
      loglik_gradS1_list, loglik_hessS1_list,
      logdiff_hess_log2plain)
    lik_hess0.list <- map3(as.list(dd_lik$M0),
      loglik_grad0_list, loglik_hess0_list,
      logdiff_hess_log2plain)
    lik_hess1.list <- map3(as.list(dd_lik$M1),
      loglik_grad1_list, loglik_hess1_list,
      logdiff_hess_log2plain)


    ### Compute, per individual, full
    ###   - Gradient vector d(loglik)/d(theta)
    ###   - Hessian matrix d2(loglik)/d(theta2)
    # Parameters for large purrr call: Likelihoods, gradients, and Hessians.
    purrr.params <- list(
      # Auxiliaries.
      as.list(dd_lik$lik), as.list(dd_lik$delta0), as.list(dd_lik$delta1),
      # Likelihoods.
      as.list(dd_lik$MS0), as.list(dd_lik$MS1),
      as.list(dd_lik$M0), as.list(dd_lik$M1),
      # Gradients.
      lik_gradS0.list, lik_gradS1.list, lik_grad0.list, lik_grad1.list,
      # Hessians.
      lik_hessS0.list, lik_hessS1.list, lik_hess0.list, lik_hess1.list)
    #  The call itself.
    loglik_grad_hess_list <- purrr::pmap(purrr.params,
      # The function that operates on all parameters.
      # Recall: Everything is valid for a singular individual.
      function(lik, delta0, delta1,
        MS0, MS1, M0, M1,
        grS0, grS1, gr0, gr1,
        heS0, heS1, he0, he1) {

        ## Likelihood gradient as a 3 x 1 block vector.
        # thetaS.
        dlik.dthetaS <- delta0 * M0 * grS0 + delta1 * M1 * grS1
        # theta0.
        dlik.dtheta0 <- delta0 * gr0 * MS0
        # theta1.
        dlik.dtheta1 <- delta1 * gr1 * MS1
        # Put together to full gradient.
        dlik.dtheta <- rbind(dlik.dthetaS, dlik.dtheta0, dlik.dtheta1)

        # Transform to log-likelihood gradient.
        loglik_grad <- logdiff_grad_plain2log(
          Fx = lik, grad_F = dlik.dtheta)
        rownames(loglik_grad) <- si_names

        ## Likelihood Hessian is a 3 x 3 block matrix.
        ## Compute only lower-triangle parts.
        # thetaS-thetaS.
        d2lik_dthetaS2 <- delta0 * M0 * heS0 + delta1 * M1 * heS1
        # theta[0|1]-thetaS.
        d2lik_dtheta0_dthetaS <- delta0 * tensor_prod(gr0, grS0)
        d2lik_dtheta1_dthetaS <- delta1 * tensor_prod(gr1, grS1)
        # theta[0|1]-theta[0|1].
        d2lik_dtheta02 <- delta0 * he0 * MS0
        d2lik_dtheta12 <- delta1 * he1 * MS1
        d2lik_dtheta1_dtheta0 <- matrix(0, nrow = d1, ncol= d0)
        # Build block matrix, one block-row at a time.
        d2lik_dtheta2_up <- cbind(d2lik_dthetaS2, t(d2lik_dtheta0_dthetaS),
          t(d2lik_dtheta1_dthetaS))
        d2lik_dtheta2_mi <- cbind(d2lik_dtheta0_dthetaS, d2lik_dtheta02,
          t(d2lik_dtheta1_dtheta0))
        d2lik_dtheta2_lo <- cbind(d2lik_dtheta1_dthetaS, d2lik_dtheta1_dtheta0,
          d2lik_dtheta12)
        # Put together to full Hessian.
        d2lik_dtheta2 <- rbind(d2lik_dtheta2_up, d2lik_dtheta2_mi, d2lik_dtheta2_lo)

        # Transform to log-likelhood Hessian.
        loglik_hess <- logdiff_hess_plain2log_hybrid(
          Fx = lik, grad_logF = loglik_grad, hess_F = d2lik_dtheta2)
        rownames(loglik_hess) <- si_names
        colnames(loglik_hess) <- si_names

        # Return both gradient (just for fun) and Hessian (needed).
        ret <- list(loglik.grad = loglik_grad, loglik.hess = loglik_hess)
        return (ret)
      }) # End of long pmap call.

    # Extract gradients and Hessians as separate lists.
    loglik_grad_list <- purrr::map(loglik_grad_hess_list, "loglik.grad")
    loglik_hess_list <- purrr::map(loglik_grad_hess_list, "loglik.hess")

    # Sum Hessians over individuals.
    loglik_hess_em <- Reduce("+", loglik_hess_list)

  } else {
    ### EM not requested, use numerical methods instead.

    ## Log-likelihood function.
    loglik_fcn <- function(theta.flat) {
      theta <- unflatten.theta(theta.flat)
      dd.lik <- setup.likelihood.data(data = data, theta = theta)
      ret <- sum(dd.lik$loglik) # Sum over individuals.
      return (ret)
    }

    ## Numerical maximization.
    res_opt <- suppressWarnings(optim(par = theta_start_flat,
      fn = loglik_fcn,
      gr = NULL, # no gradient function.
      method = "BFGS",
      hessian = F, # More accurate if computed separately.
      control=list(
        fnscale = -1, # Maximization rather than minimization.
        ndeps = rep(1e-5, length(theta_start_flat)), # Gradient step size.
        abstol = 1e-10,
        reltol = 1e-10
      )
    ))
    theta_num_flat <- res_opt$par
    theta_num <- unflatten.theta(theta_num_flat)

    # Numerical gradient (should be zero). Only for sanity checks.
    if (FALSE) {
      loglik_grad_num <- tryCatch(numDeriv::jacobian(loglik_fcn, theta_em_flat),
        warning = null.fcn, error = null.fcn)
      if (!is.null(loglik_grad_num)) {
        loglik_grad_num %<>% t # Column vector.
        rownames(loglik_grad_num) <- si_names
        max_grad_num <- max(abs(loglik_grad_num))
        stopifnot(max_grad_num < 5e-2)
      } else {
        warning("Problems with num-gradient.")
      }
    }

    # Numerical Hessian at maximum.
    loglik_hess_num <- tryCatch(numDeriv::hessian(loglik_fcn, theta_num_flat),
      warning = null.fcn, error = null.fcn)
    if (!is.null(loglik_hess_num)) {
      rownames(loglik_hess_num) <- si_names
      colnames(loglik_hess_num) <- si_names
    } else {
      stop("Problems with num-Hessian.")
    }
  } # End if (do_em)


  # Return value.
  if (do_em) {
    ret <- list(
      theta = theta_em_flat, # Estimated parameters
      info = -loglik_hess_em # Observed Fisher information.
    )
  } else {
    ret <- list(
      theta = theta_num_flat, # Estimated parameters.
      info = -loglik_hess_num # Observed Fisher information.
    )
  }

  return (ret)
} # End function latent.
