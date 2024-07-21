#' Likelihood Function for Normal Outcome Mechanism with a Binary Mediator and an Interaction Term
#'
#' @param param_start A numeric vector or column matrix of starting values for the \eqn{\theta}
#'   parameters in the outcome mechanism and \eqn{\sigma} parameter.
#'   The number of elements in \code{param_start}
#'   should be equal to the number of columns of \code{x_matrix} and \code{c_matrix} plus 2 
#'   (if \code{interaction_indicator} is \code{FALSE}) or 3 (if
#'   \code{interaction_indicator} is \code{TRUE}). Starting values should be
#'   provided in the following order: intercept, slope coefficient for the \code{x_matrix} term,
#'   slope coefficient for the mediator \code{m} term,
#'   slope coefficient for first column of the \code{c_matrix}, ...,
#'   slope coefficient for the final column of the \code{c_matrix},
#'   and, optionally, slope coefficient for \code{xm}). The final entry should be
#'   the starting value for \eqn{\sigma}.
#' @param m  vector or column matrix containing the true binary mediator or the
#'   E-step weight (with values between 0 and 1). There
#'   should be no \code{NA} terms.
#' @param x A vector or column matrix of the predictor or exposure of interest. There
#'   should be no \code{NA} terms.
#' @param c_matrix A numeric matrix of covariates in the true mediator and outcome mechanisms.
#'   \code{c_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param outcome A vector containing the outcome variables of interest. There
#'   should be no \code{NA} terms.
#' @param sample_size An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{X} or \code{Z}.
#' @param n_cat The number of categorical values that the true outcome, \code{M},
#'   and the observed outcome, \code{M*} can take.
#'
#' @return \code{theta_optim_XM} returns a numeric value of the (negative) log-likelihood function.
#'
theta_optim_XM <- function(param_start, m, x, c_matrix, outcome,
                           sample_size, n_cat){
  
  theta_v <- param_start[1:(length(param_start) - 1)]
  sigma_v <- param_start[length(param_start)]
  
  term1 <- -log(sqrt(2 * pi * sigma_v))
  
  xc_matrix <- matrix(c(rep(1, sample_size), x, c(c_matrix)),
                      nrow = sample_size, byrow = FALSE)
  theta_xc <- theta_v[-c(3, length(theta_v))]
  theta_xc_multiplied <- xc_matrix %*% theta_xc
  
  multiplicative_term <- -1 / (2 * sigma_v)
  theta_term1 <- (outcome - theta_xc_multiplied)^2
  theta_term2 <- -2 * theta_v[3] * m *
    (outcome - theta_xc_multiplied)
  theta_term3 <- theta_v[3]^2 * m
  theta_term4 <- -2 * theta_v[length(theta_v)] * x * m * (outcome - theta_xc_multiplied)
  theta_term5 <- theta_v[length(theta_v)]^2 * x^2 * m
  theta_term6 <- 2 * theta_v[3] * theta_v[length(theta_v)] * x * m
  
  summand = term1 + multiplicative_term * (theta_term1 + theta_term2 + theta_term3 +
                                             theta_term4 + theta_term5 + theta_term6)
  result = -sum(summand)
  return(result)
}
