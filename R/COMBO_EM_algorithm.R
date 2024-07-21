#' EM-Algorithm Estimation of the Binary Outcome Misclassification Model
#'
#' Jointly estimate \eqn{\beta} and \eqn{\gamma} parameters from the true outcome
#' and observation mechanisms, respectively, in a binary outcome misclassification
#' model.
#'
#' @param Ystar A numeric vector of indicator variables (1, 2) for the observed
#'   outcome \code{Y*}. There should be no \code{NA} terms. The reference category is 2.
#' @param x_matrix A numeric matrix of covariates in the true outcome mechanism.
#'   \code{x_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param z_matrix A numeric matrix of covariates in the observation mechanism.
#'   \code{z_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param beta_start A numeric vector or column matrix of starting values for the \eqn{\beta}
#'   parameters in the true outcome mechanism. The number of elements in \code{beta_start}
#'   should be equal to the number of columns of \code{x_matrix} plus 1.
#' @param gamma_start A numeric vector or matrix of starting values for the \eqn{\gamma}
#'   parameters in the observation mechanism. In matrix form, the \code{gamma_start} matrix rows
#'   correspond to parameters for the \code{Y* = 1}
#'   observed outcome, with the dimensions of \code{z_matrix} plus 1, and the
#'   gamma parameter matrix columns correspond to the true outcome categories
#'   \eqn{M \in \{1, 2\}}. A numeric vector for \code{gamma_start} is
#'   obtained by concatenating the gamma matrix, i.e. \code{gamma_start <- c(gamma_matrix)}.
#' @param tolerance A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param max_em_iterations An integer specifying the maximum number of
#'   iterations of the EM algorithm. The default is \code{1500}.
#' @param em_method A character string specifying which EM algorithm will be applied.
#'   Options are \code{"em"}, \code{"squarem"}, or \code{"pem"}. The default and
#'   recommended option is \code{"squarem"}.
#'
#' @return \code{COMBO_EM_algorithm} returns a data frame containing four columns. The first
#'   column, \code{Parameter}, represents a unique parameter value for each row.
#'   The next column contains the parameter \code{Estimates}, followed by the standard
#'   error estimates, \code{SE}. The final column, \code{Convergence}, reports
#'   whether or not the algorithm converged for a given parameter estimate.
#'
#' 
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include COMBO_weight.R
#' @include COMBO_EM_function.R
#'
#' @importFrom stats rnorm rgamma rmultinom coefficients binomial glm
#' @importFrom turboEM turboem
#' @importFrom Matrix nearPD
#'
COMBO_EM_algorithm <- function(Ystar,
                               x_matrix, z_matrix,
                               beta_start, gamma_start,
                               tolerance = 1e-7, max_em_iterations = 1500,
                               em_method = "squarem"){

  if (is.data.frame(z_matrix))
    z_matrix <- as.matrix(z_matrix)
  if (!is.numeric(z_matrix))
    stop("'z_matrix' should be a numeric matrix.")
  
  if (is.vector(z_matrix))
    z_matrix <- as.matrix(z_matrix)
  if (!is.matrix(z_matrix))
    stop("'z_matrix' should be a matrix or data.frame.")
  
  if (!is.null(x_matrix)) {
    if (is.data.frame(x_matrix))
      x_matrix <- as.matrix(x_matrix)
    if (!is.numeric(x_matrix))
      stop("'x_matrix' must be numeric.")
    if (is.vector(x_matrix))
      x_matrix <- as.matrix(x_matrix)
    if (!is.matrix(x_matrix))
      stop("'x_matrix' must be a data.frame or matrix.")
  }
  
  if (!is.numeric(Ystar) || !is.vector(Ystar))
    stop("'Ystar' must be a numeric vector.")
  if (length(setdiff(1:2, unique(Ystar))) != 0)
    stop("'Ystar' must be coded 1/2, where the reference category is 2.")
  
  n_cat = 2
  sample_size = length(Ystar)
  
  if (nrow(z_matrix) != sample_size)
    stop("The number of rows of 'z_matrix' must match the length of 'Ystar'.")
  if (!is.null(x_matrix) && nrow(x_matrix) != sample_size)
    stop("The number of rows of 'x_matrix' must match the length of 'Ystar'.")
  
  X = matrix(c(rep(1, sample_size), c(x_matrix)),
             byrow = FALSE, nrow = sample_size)
  Z = matrix(c(rep(1, sample_size), c(z_matrix)),
             byrow = FALSE, nrow = sample_size)
  
  obs_Y_reps = matrix(rep(Ystar, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix = matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                           byrow = FALSE)
  obs_Y_matrix = 1 * (obs_Y_reps == category_matrix)
  
  control_settings = list(convtype = "parameter", tol = tolerance,
                          stoptype = "maxiter", maxiter = max_em_iterations)
  
  results = turboEM::turboem(par = c(c(beta_start), c(gamma_start)),
                             fixptfn = COMBO_EM_function, 
                             method = c(em_method),
                             obs_Y_matrix = obs_Y_matrix,
                             X = X, Z = Z,
                             sample_size = sample_size, n_cat = n_cat,
                             control.run = control_settings)
  
  Ystar01 = ifelse(Ystar == 1, 1, ifelse(Ystar == 2, 0, NA))
  
  # Do label switching correction within the EM algorithm simulation
  results_i_gamma <- matrix(turboEM::pars(results)[(ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))],
                            ncol = n_cat, byrow = FALSE)
  results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)
  
  pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
  pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])
  
  flip_pistar11 <- 1 - pistar_22
  flip_pistar22 <- 1 - pistar_11
  
  J <- pistar_11 + pistar_22 - 1
  J_flip <- flip_pistar11 + flip_pistar22 - 1
  
  estimates_i <-  if ((J_flip <= J) |
                      (is.na(pistar_11) & is.na(pistar_22))) {
    # If turboem cannot estimate the parameters they will be NA.
    turboEM::pars(results)
  } else {
    gamma_index = (ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))
    n_gamma_param = length(gamma_index) / n_cat
    gamma_flip_index = ncol(X) + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)
    c(-1*turboEM::pars(results)[1:ncol(X)], turboEM::pars(results)[gamma_flip_index])
  }
  
  # Set parameter names
  beta_param_names <- paste0(rep("beta", ncol(X)), 1:ncol(X))
  gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                              rep(1:ncol(Z), n_cat),
                              rep(1:n_cat, each = ncol(Z)))
  
  estimates <- data.frame(Parameter = c(beta_param_names,
                                        gamma_param_names),
                          Estimates = c(estimates_i),
                          Convergence = c(rep(results$convergence,
                                              length(c(beta_param_names,
                                                       gamma_param_names)))))

  return(estimates)
}
