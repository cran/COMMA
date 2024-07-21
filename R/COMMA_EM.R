#' EM Algorithm Estimation of the Binary Mediator Misclassification Model
#' 
#' Jointly estimate \eqn{\beta}, \eqn{\gamma}, and \eqn{\theta} parameters from 
#' the true mediator, observed mediator, and outcome mechanisms, respectively,
#' in a binary mediator misclassification model.
#'
#' @param Mstar A numeric vector of indicator variables (1, 2) for the observed
#'   mediator \code{M*}. There should be no \code{NA} terms. The reference category is 2.
#' @param outcome A vector containing the outcome variables of interest. There
#'   should be no \code{NA} terms.
#' @param outcome_distribution A character string specifying the distribution of 
#'   the outcome variable. Options are \code{"Bernoulli"}, \code{"Normal"}, or
#'   \code{"Poisson"}.
#' @param interaction_indicator A logical value indicating if an interaction between
#'   \code{x} and \code{m} should be used to generate the outcome variable, \code{y}.
#' @param x_matrix A numeric matrix of predictors in the true mediator and outcome mechanisms.
#'   \code{x_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param z_matrix A numeric matrix of covariates in the observation mechanism.
#'   \code{z_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param c_matrix A numeric matrix of covariates in the true mediator and outcome mechanisms.
#'   \code{c_matrix} should not contain an intercept and no values should be \code{NA}.
#' @param beta_start A numeric vector or column matrix of starting values for the \eqn{\beta}
#'   parameters in the true mediator mechanism. The number of elements in \code{beta_start}
#'   should be equal to the number of columns of \code{x_matrix} and \code{c_matrix} plus 1.
#'   Starting values should be provided in the following order: intercept, slope 
#'   coefficient for the \code{x_matrix} term, slope coefficient for first column
#'   of the \code{c_matrix}, ..., slope coefficient for the final column of the \code{c_matrix}.
#' @param gamma_start A numeric vector or matrix of starting values for the \eqn{\gamma}
#'   parameters in the observation mechanism. In matrix form, the \code{gamma_start} matrix rows
#'   correspond to parameters for the \code{M* = 1}
#'   observed mediator, with the dimensions of \code{z_matrix} plus 1, and the
#'   gamma parameter matrix columns correspond to the true mediator categories
#'   \eqn{M \in \{1, 2\}}. A numeric vector for \code{gamma_start} is
#'   obtained by concatenating the gamma matrix, i.e. \code{gamma_start <- c(gamma_matrix)}.
#'   Starting values should be provided in the following order within each column:
#'   intercept, slope coefficient for first column of the \code{z_matrix}, ...,
#'   slope coefficient for the final column of the \code{z_matrix}.
#' @param theta_start A numeric vector or column matrix of starting values for the \eqn{\theta}
#'   parameters in the outcome mechanism. The number of elements in \code{theta_start}
#'   should be equal to the number of columns of \code{x_matrix} and \code{c_matrix} plus 2 
#'   (if \code{interaction_indicator} is \code{FALSE}) or 3 (if
#'   \code{interaction_indicator} is \code{TRUE}). Starting values should be
#'   provided in the following order: intercept, slope coefficient for the \code{x_matrix} term,
#'   slope coefficient for the mediator \code{m} term,
#'   slope coefficient for first column of the \code{c_matrix}, ...,
#'   slope coefficient for the final column of the \code{c_matrix},
#'   and, optionally, slope coefficient for \code{xm}).
#' @param sigma_start A numeric value specifying the starting value for the
#'   standard deviation. This value is only required if \code{outcome_distribution}
#'   is \code{"Normal"}. Otherwise, this value is set to \code{NULL}.
#' @param tolerance A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param max_em_iterations A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param em_method A character string specifying which EM algorithm will be applied.
#'   Options are \code{"em"}, \code{"squarem"}, or \code{"pem"}. The default and
#'   recommended option is \code{"squarem"}.
#'
#' @return \code{COMMA_EM} returns a data frame containing four columns. The first
#'   column, \code{Parameter}, represents a unique parameter value for each row.
#'   The next column contains the parameter \code{Estimates}, followed by the standard
#'   error estimates, \code{SE}. The final column, \code{Convergence}, reports
#'   whether or not the algorithm converged for a given parameter estimate.
#'   
#' @export
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include w_m_binaryY.R
#' @include w_m_normalY.R
#' @include EM_function_bernoulliY.R
#' @include EM_function_bernoulliY_XM.R
#' @include EM_function_normalY.R
#' @include EM_function_normalY_XM.R
#' @include COMMA_data.R
#'
#' @importFrom stats rnorm rgamma rmultinom coefficients binomial glm
#' @importFrom turboEM turboem
#' @importFrom Matrix nearPD
#'
#' @examples
#' set.seed(20240709)
#' sample_size <- 2000
#' 
#' n_cat <- 2 # Number of categories in the binary mediator
#' 
#' # Data generation settings
#' x_mu <- 0
#' x_sigma <- 1
#' z_shape <- 1
#' c_shape <- 1
#' 
#' # True parameter values (gamma terms set the misclassification rate)
#' true_beta <- matrix(c(1, -2, .5), ncol = 1)
#' true_gamma <- matrix(c(1, 1, -.5, -1.5), nrow = 2, byrow = FALSE)
#' true_theta <- matrix(c(1, 1.5, -2, -.2), ncol = 1)
#' 
#' example_data <- COMMA_data(sample_size, x_mu, x_sigma, z_shape, c_shape,
#'                            interaction_indicator = FALSE,
#'                            outcome_distribution = "Bernoulli",
#'                            true_beta, true_gamma, true_theta)
#'                            
#' beta_start <- matrix(rep(1, 3), ncol = 1)
#' gamma_start <- matrix(rep(1, 4), nrow = 2, ncol = 2)
#' theta_start <- matrix(rep(1, 4), ncol = 1)
#' 
#' Mstar = example_data[["obs_mediator"]]
#' outcome = example_data[["outcome"]]
#' x_matrix = example_data[["x"]]
#' z_matrix = example_data[["z"]]
#' c_matrix = example_data[["c"]]
#'                            
#' EM_results <- COMMA_EM(Mstar, outcome, "Bernoulli", FALSE,
#'                        x_matrix, z_matrix, c_matrix,
#'                        beta_start, gamma_start, theta_start)
#' 
#' EM_results
#'
COMMA_EM <- function(Mstar, # Observed mediator vector
                     outcome, # Outcome vector
                     outcome_distribution,
                     interaction_indicator,
                     # Predictor matrices
                     x_matrix, z_matrix, c_matrix,
                     # Start values for parameters
                     beta_start, gamma_start,
                     theta_start, sigma_start = NULL,
                     # EM settings
                     tolerance = 1e-7, max_em_iterations = 1500,
                     em_method = "squarem"){

  n_cat = 2 # Number of categories in mediator
  sample_size = length(Mstar) # Sample size
  
  # Create design matrices
  X = matrix(c(rep(1, sample_size), c(x_matrix)),
             byrow = FALSE, nrow = sample_size)
  Z = matrix(c(rep(1, sample_size), c(z_matrix)),
             byrow = FALSE, nrow = sample_size)
  
  x_mat <- as.matrix(x_matrix)
  c_mat <- as.matrix(c_matrix)
  
  # Create a matrix of observed mediator variables using dummy coding
  obs_M_reps = matrix(rep(Mstar, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix = matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                           byrow = FALSE)
  obs_M_matrix = 1 * (obs_M_reps == category_matrix)
  
  # EM algorithm settings
  control_settings = list(convtype = "parameter", tol = tolerance,
                          stoptype = "maxiter", maxiter = max_em_iterations)
  
  # Run EM algorithm using turboEM
  
  if(interaction_indicator == FALSE & outcome_distribution == "Bernoulli"){
    
    results = turboEM::turboem(par = c(c(beta_start), c(gamma_start),
                                       c(theta_start)),
                               fixptfn = EM_function_bernoulliY,
                               method = c(em_method),
                               obs_mediator = Mstar,
                               obs_outcome = outcome,
                               X = X, Z = Z, c_matrix = c_matrix,
                               sample_size = sample_size, n_cat = n_cat,
                               control.run = control_settings)
    
    # Recode observed mediator from 1/2, make 1/0
    Mstar01 = ifelse(Mstar == 1, 1, ifelse(Mstar == 2, 0, NA))
    
    # Do label switching correction within the EM algorithm simulation
    design_matrix = cbind(X, c_matrix)
    
    gamma_index_1 = ncol(design_matrix) + 1
    gamma_index_2 = gamma_index_1 + (ncol(Z) * 2) - 1
    
    n_param <- length(turboEM::pars(results))
    
    results_i_gamma <- matrix(turboEM::pars(results)[gamma_index_1:gamma_index_2],
                              ncol = n_cat, byrow = FALSE)
    results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)
    
    pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
    pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])
    
    flip_pistar11 <- 1 - pistar_22
    flip_pistar22 <- 1 - pistar_11
    
    J <- pistar_11 + pistar_22 - 1
    J_flip <- flip_pistar11 + flip_pistar22 - 1
    
    estimates_i <- if ((J_flip <= J) |
                       (is.na(pistar_11) & is.na(pistar_22))) {
      # If turboem cannot estimate the parameters they will be NA.
      turboEM::pars(results)
    } else {
      gamma_index = gamma_index_1:gamma_index_2
      n_gamma_param = length(gamma_index) / n_cat
      gamma_flip_index = ncol(design_matrix) + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)
      
      c(-1*turboEM::pars(results)[1:ncol(design_matrix)], # beta * -1
        turboEM::pars(results)[gamma_flip_index], # flip gammas
        turboEM::pars(results)[gamma_index_2 + 1] + turboEM::pars(results)[n_param - ncol(c_mat)], # add theta_m to intercept
        turboEM::pars(results)[(gamma_index_2 + 2):(gamma_index_2 + 1 + ncol(x_mat))],
        -1 * turboEM::pars(results)[gamma_index_2 + 1 + ncol(x_mat) + 1], # multiply theta_m by -1 
        turboEM::pars(results)[(gamma_index_2 + 1 + ncol(x_mat) + 1 + 1):n_param])
    }
    
    # Set parameter names
    beta_param_names <- paste0(rep("beta_", ncol(design_matrix)), 0:(ncol(design_matrix) - 1))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           paste0(rep("theta_x", ncol(x_mat)), 1:ncol(x_mat)),
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)))
    
    # Return data frame of estimates
    estimates <- data.frame(Parameter = c(beta_param_names,
                                          gamma_param_names,
                                          theta_param_names),
                            Estimates = c(estimates_i),
                            Convergence = c(rep(results$convergence,
                                                length(c(beta_param_names,
                                                         gamma_param_names,
                                                         theta_param_names)))))
    
  } else if(interaction_indicator == TRUE & outcome_distribution == "Bernoulli"){
    
    results = turboEM::turboem(par = c(c(beta_start), c(gamma_start),
                                       c(theta_start)),
                               fixptfn = EM_function_bernoulliY_XM,
                               method = c(em_method),
                               obs_mediator = Mstar,
                               obs_outcome = outcome,
                               X = X, Z = Z, c_matrix = c_matrix,
                               sample_size = sample_size, n_cat = n_cat,
                               control.run = control_settings)
    
    # Recode observed mediator from 1/2, make 1/0
    Mstar01 = ifelse(Mstar == 1, 1, ifelse(Mstar == 2, 0, NA))
    
    # Do label switching correction within the EM algorithm simulation
    design_matrix = cbind(X, c_matrix)
    
    gamma_index_1 = ncol(design_matrix) + 1
    gamma_index_2 = gamma_index_1 + (ncol(Z) * 2) - 1
    
    n_param <- length(turboEM::pars(results))
    
    results_i_gamma <- matrix(turboEM::pars(results)[gamma_index_1:gamma_index_2],
                              ncol = n_cat, byrow = FALSE)
    results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)
    
    pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
    pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])
    
    flip_pistar11 <- 1 - pistar_22
    flip_pistar22 <- 1 - pistar_11
    
    J <- pistar_11 + pistar_22 - 1
    J_flip <- flip_pistar11 + flip_pistar22 - 1
    
    estimates_i <- if ((J_flip <= J) |
                       (is.na(pistar_11) & is.na(pistar_22))) {
      # If turboem cannot estimate the parameters they will be NA.
      turboEM::pars(results)
    } else {
      gamma_index = gamma_index_1:gamma_index_2
      n_gamma_param = length(gamma_index) / n_cat
      gamma_flip_index = ncol(design_matrix) + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)
      
      c(-1*turboEM::pars(results)[1:ncol(design_matrix)], # beta * -1
        turboEM::pars(results)[gamma_flip_index], # flip gammas
        turboEM::pars(results)[gamma_index_2 + 1] + turboEM::pars(results)[n_param - ncol(c_mat)], # add theta_m to intercept
        turboEM::pars(results)[(gamma_index_2 + 2):(gamma_index_2 + 1 + ncol(x_mat))] + turboEM::pars(results)[n_param],
        -1 * turboEM::pars(results)[gamma_index_2 + 1 + ncol(x_mat) + 1], # multiply theta_m by -1 
        turboEM::pars(results)[(gamma_index_2 + 1 + ncol(x_mat) + 1 + 1):(n_param - 1)],
        -1 * turboEM::pars(results)[n_param])
    }
    
    # Set parameter names
    beta_param_names <- paste0(rep("beta_", ncol(design_matrix)), 0:(ncol(design_matrix) - 1))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           "theta_x",
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)),
                           "theta_xm")
    
    # Return data frame of estimates
    estimates <- data.frame(Parameter = c(beta_param_names,
                                          gamma_param_names,
                                          theta_param_names),
                            Estimates = c(estimates_i),
                            Convergence = c(rep(results$convergence,
                                                length(c(beta_param_names,
                                                         gamma_param_names,
                                                         theta_param_names)))))
    
  } else if(interaction_indicator == FALSE & outcome_distribution == "Normal"){
    
    results = turboEM::turboem(par = c(c(beta_start), c(gamma_start),
                                       c(theta_start), sigma_start),
                               fixptfn = EM_function_normalY,
                               method = c(em_method),
                               obs_mediator = Mstar,
                               obs_outcome = outcome,
                               X = X, Z = Z, c_matrix = c_matrix,
                               sample_size = sample_size, n_cat = n_cat,
                               control.run = control_settings)
    
    # Recode observed mediator from 1/2, make 1/0
    Mstar01 = ifelse(Mstar == 1, 1, ifelse(Mstar == 2, 0, NA))
    
    # Do label switching correction within the EM algorithm simulation
    design_matrix = cbind(X, c_matrix)
    
    gamma_index_1 = ncol(design_matrix) + 1
    gamma_index_2 = gamma_index_1 + (ncol(Z) * 2) - 1
    
    n_param <- length(turboEM::pars(results))
    
    results_i_gamma <- matrix(turboEM::pars(results)[gamma_index_1:gamma_index_2],
                              ncol = n_cat, byrow = FALSE)
    results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)
    
    pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
    pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])
    
    flip_pistar11 <- 1 - pistar_22
    flip_pistar22 <- 1 - pistar_11
    
    J <- pistar_11 + pistar_22 - 1
    J_flip <- flip_pistar11 + flip_pistar22 - 1
    
    estimates_i <- if ((J_flip <= J) |
                       (is.na(pistar_11) & is.na(pistar_22))) {
      # If turboem cannot estimate the parameters they will be NA.
      turboEM::pars(results)
    } else {
      gamma_index = gamma_index_1:gamma_index_2
      n_gamma_param = length(gamma_index) / n_cat
      gamma_flip_index = ncol(design_matrix) + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)
      
      c(-1*turboEM::pars(results)[1:ncol(design_matrix)], # beta * -1
        turboEM::pars(results)[gamma_flip_index], # flip gammas
        turboEM::pars(results)[gamma_index_2 + 1] + turboEM::pars(results)[n_param - ncol(c_mat)], # add theta_m to intercept
        turboEM::pars(results)[(gamma_index_2 + 2):(gamma_index_2 + 1 + ncol(x_mat))],
        -1 * turboEM::pars(results)[gamma_index_2 + 1 + ncol(x_mat) + 1], # multiply theta_m by -1 
        turboEM::pars(results)[(gamma_index_2 + 1 + ncol(x_mat) + 1 + 1):n_param])
    }
    
    # Set parameter names
    beta_param_names <- paste0(rep("beta_", ncol(design_matrix)), 0:(ncol(design_matrix) - 1))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           paste0(rep("theta_x", ncol(x_mat)), 1:ncol(x_mat)),
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)))
    
    # Return data frame of estimates
    estimates <- data.frame(Parameter = c(beta_param_names,
                                          gamma_param_names,
                                          theta_param_names,
                                          "sigma"),
                            Estimates = c(estimates_i),
                            Convergence = c(rep(results$convergence,
                                                length(c(beta_param_names,
                                                         gamma_param_names,
                                                         theta_param_names,
                                                         "sigma")))))
    
  } else if(interaction_indicator == TRUE & outcome_distribution == "Normal"){
    
    results = turboEM::turboem(par = c(c(beta_start), c(gamma_start),
                                       c(theta_start), sigma_start),
                               fixptfn = EM_function_normalY_XM,
                               method = c(em_method),
                               obs_mediator = Mstar,
                               obs_outcome = outcome,
                               X = X, Z = Z, c_matrix = c_matrix,
                               sample_size = sample_size, n_cat = n_cat,
                               control.run = control_settings)
    
    # Recode observed mediator from 1/2, make 1/0
    Mstar01 = ifelse(Mstar == 1, 1, ifelse(Mstar == 2, 0, NA))
    
    # Do label switching correction within the EM algorithm simulation
    design_matrix = cbind(X, c_matrix)
    
    gamma_index_1 = ncol(design_matrix) + 1
    gamma_index_2 = gamma_index_1 + (ncol(Z) * 2) - 1
    
    n_param <- length(turboEM::pars(results))
    
    results_i_gamma <- matrix(turboEM::pars(results)[gamma_index_1:gamma_index_2],
                              ncol = n_cat, byrow = FALSE)
    results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)
    
    pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
    pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])
    
    flip_pistar11 <- 1 - pistar_22
    flip_pistar22 <- 1 - pistar_11
    
    J <- pistar_11 + pistar_22 - 1
    J_flip <- flip_pistar11 + flip_pistar22 - 1
    
    estimates_i <- if ((J_flip <= J) |
                       (is.na(pistar_11) & is.na(pistar_22))) {
      # If turboem cannot estimate the parameters they will be NA.
      turboEM::pars(results)
    } else {
      gamma_index = gamma_index_1:gamma_index_2
      n_gamma_param = length(gamma_index) / n_cat
      gamma_flip_index = ncol(design_matrix) + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)
      
      c(-1*turboEM::pars(results)[1:ncol(design_matrix)], # beta * -1
        turboEM::pars(results)[gamma_flip_index], # flip gammas
        turboEM::pars(results)[gamma_index_2 + 1] + turboEM::pars(results)[n_param - ncol(c_mat)], # add theta_m to intercept
        turboEM::pars(results)[(gamma_index_2 + 2):(gamma_index_2 + 1 + ncol(x_mat))] + turboEM::pars(results)[n_param],
        -1 * turboEM::pars(results)[gamma_index_2 + 1 + ncol(x_mat) + 1], # multiply theta_m by -1 
        turboEM::pars(results)[(gamma_index_2 + 1 + ncol(x_mat) + 1 + 1):(n_param - 2)],
        -1 * turboEM::pars(results)[n_param - 1], # multiply interaction term by -1
        turboEM::pars(results)[n_param]) # Make sure you don't multiply the sigma!
    }
    
    # Set parameter names
    beta_param_names <- paste0(rep("beta_", ncol(design_matrix)), 0:(ncol(design_matrix) - 1))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           "theta_x",
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)),
                           "theta_xm")
    
    # Return data frame of estimates
    estimates <- data.frame(Parameter = c(beta_param_names,
                                          gamma_param_names,
                                          theta_param_names,
                                          "sigma"),
                            Estimates = c(estimates_i),
                            Convergence = c(rep(results$convergence,
                                                length(c(beta_param_names,
                                                         gamma_param_names,
                                                         theta_param_names,
                                                         "sigma")))))
    
  } else if(interaction_indicator == FALSE & outcome_distribution == "Poisson"){
    
    results = turboEM::turboem(par = c(c(beta_start), c(gamma_start),
                                       c(theta_start)),
                               fixptfn = EM_function_poissonY,
                               method = c(em_method),
                               obs_mediator = Mstar,
                               obs_outcome = outcome,
                               X = X, Z = Z, c_matrix = c_matrix,
                               sample_size = sample_size, n_cat = n_cat,
                               control.run = control_settings)
    
    # Recode observed mediator from 1/2, make 1/0
    Mstar01 = ifelse(Mstar == 1, 1, ifelse(Mstar == 2, 0, NA))
    
    # Do label switching correction within the EM algorithm simulation
    design_matrix = cbind(X, c_matrix)
    
    gamma_index_1 = ncol(design_matrix) + 1
    gamma_index_2 = gamma_index_1 + (ncol(Z) * 2) - 1
    
    n_param <- length(turboEM::pars(results))
    
    results_i_gamma <- matrix(turboEM::pars(results)[gamma_index_1:gamma_index_2],
                              ncol = n_cat, byrow = FALSE)
    results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)
    
    pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
    pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])
    
    flip_pistar11 <- 1 - pistar_22
    flip_pistar22 <- 1 - pistar_11
    
    J <- pistar_11 + pistar_22 - 1
    J_flip <- flip_pistar11 + flip_pistar22 - 1
    
    estimates_i <- if ((J_flip <= J) |
                       (is.na(pistar_11) & is.na(pistar_22))) {
      # If turboem cannot estimate the parameters they will be NA.
      turboEM::pars(results)
    } else {
      gamma_index = gamma_index_1:gamma_index_2
      n_gamma_param = length(gamma_index) / n_cat
      gamma_flip_index = ncol(design_matrix) + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)
      
      c(-1*turboEM::pars(results)[1:ncol(design_matrix)], # beta * -1
        turboEM::pars(results)[gamma_flip_index], # flip gammas
        turboEM::pars(results)[gamma_index_2 + 1] + turboEM::pars(results)[n_param - ncol(c_mat)], # add theta_m to intercept
        turboEM::pars(results)[(gamma_index_2 + 2):(gamma_index_2 + 1 + ncol(x_mat))],
        -1 * turboEM::pars(results)[gamma_index_2 + 1 + ncol(x_mat) + 1], # multiply theta_m by -1 
        turboEM::pars(results)[(gamma_index_2 + 1 + ncol(x_mat) + 1 + 1):n_param])
    }
    
    # Set parameter names
    beta_param_names <- paste0(rep("beta_", ncol(design_matrix)), 0:(ncol(design_matrix) - 1))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           paste0(rep("theta_x", ncol(x_mat)), 1:ncol(x_mat)),
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)))
    
    # Return data frame of estimates
    estimates <- data.frame(Parameter = c(beta_param_names,
                                          gamma_param_names,
                                          theta_param_names),
                            Estimates = c(estimates_i),
                            Convergence = c(rep(results$convergence,
                                                length(c(beta_param_names,
                                                         gamma_param_names,
                                                         theta_param_names)))))
    
  } else if(interaction_indicator == TRUE & outcome_distribution == "Poisson"){
    
    results = turboEM::turboem(par = c(c(beta_start), c(gamma_start),
                                       c(theta_start)),
                               fixptfn = EM_function_poissonY_XM,
                               method = c(em_method),
                               obs_mediator = Mstar,
                               obs_outcome = outcome,
                               X = X, Z = Z, c_matrix = c_matrix,
                               sample_size = sample_size, n_cat = n_cat,
                               control.run = control_settings)
    
    # Recode observed mediator from 1/2, make 1/0
    Mstar01 = ifelse(Mstar == 1, 1, ifelse(Mstar == 2, 0, NA))
    
    # Do label switching correction within the EM algorithm simulation
    design_matrix = cbind(X, c_matrix)
    
    gamma_index_1 = ncol(design_matrix) + 1
    gamma_index_2 = gamma_index_1 + (ncol(Z) * 2) - 1
    
    n_param <- length(turboEM::pars(results))
    
    results_i_gamma <- matrix(turboEM::pars(results)[gamma_index_1:gamma_index_2],
                              ncol = n_cat, byrow = FALSE)
    results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)
    
    pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
    pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])
    
    flip_pistar11 <- 1 - pistar_22
    flip_pistar22 <- 1 - pistar_11
    
    J <- pistar_11 + pistar_22 - 1
    J_flip <- flip_pistar11 + flip_pistar22 - 1
    
    estimates_i <- if ((J_flip <= J) |
                       (is.na(pistar_11) & is.na(pistar_22))) {
      # If turboem cannot estimate the parameters they will be NA.
      turboEM::pars(results)
    } else {
      gamma_index = gamma_index_1:gamma_index_2
      n_gamma_param = length(gamma_index) / n_cat
      gamma_flip_index = ncol(design_matrix) + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)
      
      c(-1*turboEM::pars(results)[1:ncol(design_matrix)], # beta * -1
        turboEM::pars(results)[gamma_flip_index], # flip gammas
        turboEM::pars(results)[gamma_index_2 + 1] + turboEM::pars(results)[n_param - ncol(c_mat)], # add theta_m to intercept
        turboEM::pars(results)[(gamma_index_2 + 2):(gamma_index_2 + 1 + ncol(x_mat))] + turboEM::pars(results)[n_param],
        -1 * turboEM::pars(results)[gamma_index_2 + 1 + ncol(x_mat) + 1], # multiply theta_m by -1 
        turboEM::pars(results)[(gamma_index_2 + 1 + ncol(x_mat) + 1 + 1):(n_param - 1)],
        -1 * turboEM::pars(results)[n_param])
    }
    
    # Set parameter names
    beta_param_names <- paste0(rep("beta_", ncol(design_matrix)), 0:(ncol(design_matrix) - 1))
    gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                                rep(1:ncol(Z), n_cat),
                                rep(1:n_cat, each = ncol(Z)))
    theta_param_names <- c("theta_0",
                           "theta_x",
                           "theta_m",
                           paste0(rep("theta_c", ncol(c_mat)), 1:ncol(c_mat)),
                           "theta_xm")
    
    # Return data frame of estimates
    estimates <- data.frame(Parameter = c(beta_param_names,
                                          gamma_param_names,
                                          theta_param_names),
                            Estimates = c(estimates_i),
                            Convergence = c(rep(results$convergence,
                                                length(c(beta_param_names,
                                                         gamma_param_names,
                                                         theta_param_names)))))
    
  } else {
    "Undefined"
  }
  

  return(estimates)
}
