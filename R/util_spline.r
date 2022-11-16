# Extracts inner knots.
spline_inner <- function(k) k[-c(1, length(k))]
# Extracts boundary knots.
spline_boundary <- function(k) k[c(1, length(k))]

# Cubic basis function [(x - a)_+]^3 and its derivative.
p3 <- function(x, a, deriv) {
  if (deriv) {
    ret <- 3 * pmax(x - a, 0)^2
  } else {
    ret <- pmax(x - a, 0)^3
  }
  return (ret)
}

# Spline basis functions bj(x), 0 <= j <= m + 1,
# where m = the number of internal knots.
#
# b0(x) = 1, b1(x) = x, b_j(x) = v_{j-1}(x) for 2 <= j <= m + 1,
# where v are the basis functions from flexsurv.
spline_basis_fcn <- function(x, j, knots, deriv) {
  stopifnot(is.element(j, seq_along(knots) - 1)) # 0 <= j <= m + 1.
  if (j == 0) {
    if (deriv) {
      ret <-  rep(0, length(x))
    } else {
      ret <-  rep(1, length(x))
    }
  } else if (j == 1) {
    if (deriv) {
      ret <-  rep(1, length(x))
    } else {
      ret <- x
    }
  } else {
    # j >= 2. Apply v_{j-1}.
    j1 <- j - 1
    k_min <- knots %>% head(1)
    k_int <- knots[-c(1, length(knots))]
    k_j1 <- k_int[j1]
    k_max <- knots %>% tail(1)
    lambda_j1 <- (k_max - k_j1) / (k_max - k_min)
    ret <- p3(x, k_j1, deriv) -
      lambda_j1 * p3(x, k_min, deriv) -
      (1 - lambda_j1) * p3(x, k_max, deriv)
 }
  return (ret)
}

# Array of spline basis functions (b_0, ..., b_{m+1}).
# Returned in the form of a matrix, with the b's as columns.
spline_basis_fcn_arr <- function(x, knots, deriv = F) {
  ret <- matrix(NA, nrow = length(x), ncol = length(knots))
  for (j in seq_along(knots) - 1) {
    ret[, j + 1] <- spline_basis_fcn(x, j, knots, deriv)
  }
  return (ret)
}

# Entire spline function = linear combination of basis function.
spline_fcn <- function(x, gamma, knots, deriv = F) {
  stopifnot(length(gamma) == length(knots))
  bv <- spline_basis_fcn_arr(x, knots, deriv)
  ret <- bv %*% cbind(gamma)
  ret <- ret[, 1] # Un-matrix.
  return (ret)
}

# Convenience wrapper.
spline_fcn_deriv <- function(...) spline_fcn(..., deriv = T)
