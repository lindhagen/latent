### Tensor stuff.
#' @title Tensor product
#' @description Forms the tensor product between two vectors.
#' @param x,y Two vectors of the same length, given either as plain arrays
#' or as column \code{matrix} objects.
#' @return The tensor product of \code{x} and \code{y} as a square
#' \code{matrix}.
#' @export
tensor_prod <- function(x, y) {
  ret <- x %*% t(y)
  return (ret)
}
#' @title Tensor square
#' @description Forms the tensor product of a vector with itself.
#' @param x A vector.
#' @return The tensor product (see \link{tensor_prod}) between x and itself.
#' @export
tensor_square <- function(x) {
  ret <- tensor_prod(x, x)
  return (ret)
}


### Logarithmic differentiation. Given a function F, these R functions
### convert between gradient/Hessian of F ("plain") and of log(F) ("log").
### Optionally, the Hessian can be skipped (set to NULL).
### In both cases, the function value itself (Fx = F(x)) is needed too.

## Derivatives of log(F) --> Derivatives of F.
# Gradient.
logdiff_grad_log2plain <- function(Fx, grad_logF) {
  grad_F <- Fx * grad_logF
  return (grad_F)
}
# Hessian.
logdiff_hess_log2plain <- function(Fx, grad_logF, hess_logF) {
  hess_F = Fx * (hess_logF + tensor_square(grad_logF))
  return (hess_F)
}
# Like above, but hybrid version, where the gradient of F it given,
# rather than that of log(F).
logdiff_hess_log2plain_hybrid <- function(Fx, grad_F, hess_logF) {
  hess_F = Fx * (hess_logF + tensor_square(grad_F / Fx))
  return (hess_F)
}

## Derivatives of F --> Derivatives of log(F).
# Gradient.
logdiff_grad_plain2log <- function(Fx, grad_F) {
  grad_logF <- grad_F / Fx
  return (grad_logF)
}
# Hessian.
logdiff_hess_plain2log <- function(Fx, grad_F, hess_F) {
  hess_logF <- hess_F / Fx - tensor_square(grad_F / Fx)
  return (hess_logF)
}
# Like above, but hybrid version, where the gradient of log(F) it given,
# rather than that of F.
logdiff_hess_plain2log_hybrid <- function(Fx, grad_logF, hess_F) {
  hess_logF <- hess_F / Fx - tensor_square(grad_logF)
  return (hess_logF)
}

# purrr(ish) functions.
map2 <- purrr::map2
map3 <- function(l1, l2, l3, f) {
  purrr::pmap(list(l1, l2, l3), f)
}


# Flattens a list of lists, returning it as a plain array with good names.
#' @title Flatten a list of lists
#' @description Converts a list of lists to a single list. Useful for parameter
#' handling.
#' @param lst A list consisting of lists. These lists should in turn consist of
#' either a single number, or a named array.
#' @return An array consisting of all the elements of the above lists.
#' @export
flatten.listlist <- function(lst) {
  ret <- c()
  for (xn in names(lst)) {
    x <- lst[[xn]]
    if (is.list(x)) {
      x <- unlist(x)
      names(x) <- sprintf("%s.%s", xn, names(x))
    } else {
      names(x) <- xn
    }
    ret <- c(ret, x)
  }
  return (ret)
}

# Inverse of above function, tailored for our pi-theta0-theta1 structure.
unflatten.theta <- function(theta.flat) {
  ixS <- grep("^thetaS", names(theta.flat))
  ix0 <- grep("^theta0", names(theta.flat))
  ix1 <- grep("^theta1", names(theta.flat))
  stopifnot(identical(c(ixS, ix0, ix1), seq_along(theta.flat)))

  thetaS <- theta.flat[ixS]
  theta0 <- theta.flat[ix0]
  theta1 <- theta.flat[ix1]
  names(thetaS) <- gsub("^thetaS\\.", "", names(thetaS))
  names(theta0) <- gsub("^theta0\\.", "", names(theta0))
  names(theta1) <- gsub("^theta1\\.", "", names(theta1))

  ret <- list(
    thetaS = as.list(thetaS),
    theta0 = as.list(theta0),
    theta1 = as.list(theta1)
  )
  return (ret)
}
