#' Ordinary Least Squares Estimation of the Binary Mediator Misclassification Model
#' 
#' Estimate \eqn{\beta}, \eqn{\gamma}, and \eqn{\theta} parameters from 
#' the true mediator, observed mediator, and outcome mechanisms, respectively,
#' in a binary mediator misclassification model using an ordinary least squares
#' correction.
#' 
#' Note that this method can only be used for Normal outcome models, and interaction
#' terms (between \code{x} and \code{m}) are not supported.
#'
#' @param Mstar A numeric vector of indicator variables (1, 2) for the observed
#'   mediator \code{M*}. There should be no \code{NA} terms. The reference category is 2.
#' @param outcome A vector containing the outcome variables of interest. There
#'   should be no \code{NA} terms.
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
#'   should be equal to the number of columns of \code{x_matrix} and \code{c_matrix}
#'   plus 2. Starting values should be provided in the following order: intercept,
#'   slope coefficient for the \code{x_matrix} term, slope coefficient for the
#'   mediator \code{m} term, slope coefficient for first column of the \code{c_matrix}, ...,
#'   slope coefficient for the final column of the \code{c_matrix}.
#' @param tolerance A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param max_em_iterations A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param em_method A character string specifying which EM algorithm will be applied.
#'   Options are \code{"em"}, \code{"squarem"}, or \code{"pem"}. The default and
#'   recommended option is \code{"squarem"}.
#'
#' @return \code{COMMA_PVW} returns a data frame containing four columns. The first
#'   column, \code{Parameter}, represents a unique parameter value for each row.
#'   The next column contains the parameter \code{Estimates}. The third column,
#'   \code{Convergence}, reports whether or not the algorithm converged for a
#'   given parameter estimate. The final column, \code{Method}, reports
#'   that the estimates are obtained from the "PVW" procedure. 
#'   
#' @export
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include COMBO_EM_algorithm.R
#' @include COMBO_EM_function.R
#' @include COMBO_weight.R
#' @include COMMA_data.R
#' 
#' @importFrom stats cov median
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
#' true_theta <- matrix(c(1, 1.5, -2, 2), ncol = 1)
#' 
#' example_data <- COMMA_data(sample_size, x_mu, x_sigma, z_shape, c_shape,
#'                            interaction_indicator = FALSE,
#'                            outcome_distribution = "Normal",
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
#' OLS_results <- COMMA_OLS(Mstar, outcome,
#'                          x_matrix, z_matrix, c_matrix,
#'                          beta_start, gamma_start, theta_start)
#'                          
#' OLS_results
#'
#' 
COMMA_OLS <- function(Mstar, # Observed mediator vector
                      outcome, # Outcome vector
                      # Predictor matrices
                      x_matrix, z_matrix, c_matrix,
                      # Start values for parameters
                      beta_start, gamma_start,
                      theta_start,
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
  
  # Create matrix of true mediation model predictors
  mediation_model_predictors <- cbind(x_matrix, c_matrix)
  
  # Run the COMBO EM algorithm for the true and observed mediation model
  COMBO_EM_results <- COMBO_EM_algorithm(Mstar,
                                         mediation_model_predictors,
                                         z_matrix,
                                         beta_start, gamma_start,
                                         tolerance, max_em_iterations,
                                         em_method)
  # Save results
  gamma_index_1 = ncol(mediation_model_predictors) + 2
  gamma_index_2 = gamma_index_1 + (ncol(Z) * 2) - 1
  
  n_param <- length(c(beta_start, c(gamma_start), theta_start))
  
  predicted_beta <- matrix(COMBO_EM_results$Estimates[1:(ncol(mediation_model_predictors) + 1)],
                           ncol = 1)
  predicted_gamma <- matrix(COMBO_EM_results$Estimates[gamma_index_1:gamma_index_2], 
                            ncol = 2, byrow = FALSE)
  
  # Create a matrix of observed mediator variables using dummy coding
  mstar_matrix <- matrix(c(ifelse(Mstar == 1, 1, 0),
                           ifelse(Mstar == 2, 1, 0)),
                         ncol = 2, byrow = FALSE)
  
  # Create matrix of predictors for the true mediator
  X_design <- cbind(rep(1, sample_size), mediation_model_predictors)
  
  # Generate probabilities for the true mediator value based on EM results
  pi_matrix <- pi_compute(predicted_beta, X_design, sample_size, n_cat)
  
  # Create matrix of predictors for the observed mediator
  Z_design <- matrix(c(rep(1, sample_size), z_matrix),
                     nrow = sample_size, byrow = FALSE)
  
  # Generate probabilities for observed mediator conditional on true mediator
  ## Based on EM results
  pistar_matrix <- pistar_compute(predicted_gamma, Z_design, sample_size, n_cat)
  
  # Estimate sensitivity and specificity
  sensitivity <- pistar_matrix[1:sample_size, 1]
  specificity <- pistar_matrix[(sample_size + 1):(2 * sample_size), 2]
  
  # Compute the observed mediator prevalence
  prevalence <- length(which(Mstar == 1)) / sample_size
  
  # Compute average misclassification rates
  pistar12 <- pistar_matrix[1:sample_size, 2]
  pistar21 <- pistar_matrix[(sample_size + 1):(2 * sample_size), 1]
  
  # Compute correction parameters from Nguimkeu, Rosenman, and Tennekoon (2021)
  theta_Nguimkeu <- (pistar12 + pistar21) / (1 - pistar12 - pistar21)
  squiggle_Nguimkeu <- 1 - (((prevalence - pistar12)*(1 - pistar21 - prevalence)) / 
                              ((1 - pistar12 - pistar21)*(1 - prevalence)*prevalence))
  
  # Compute covariances for the correction
  m_matrix <- matrix(ifelse(Mstar == 1, 1, 0), ncol = 1)
  sd_dd <- cov(m_matrix)
  
  predictor_matrix <- matrix(c(x_matrix, c_matrix), nrow = sample_size, byrow = FALSE)
  sd_xx <- cov(predictor_matrix)
  
  sd_xd <- cov(predictor_matrix, m_matrix)
  sd_dx <- cov(m_matrix, predictor_matrix)
  
  y_matrix <- matrix(outcome, ncol = 1)
  sd_yd <- cov(y_matrix, m_matrix)
  sd_yx <- cov(y_matrix, predictor_matrix)
  
  n_theta <- length(theta_start)
  block1_dd <- (1 - median(squiggle_Nguimkeu)) * sd_dd[1,1]
  block1_xd <- (1 + median(theta_Nguimkeu)) * sd_xd
  
  block_1_matrix <- matrix(rep(NA, (n_theta - 1)^2),
                           nrow = n_theta - 1, ncol = n_theta - 1)
  block_1_matrix[,1] <- c(block1_dd, block1_xd)
  for(i in 2:(n_theta - 1)){
    block_1_matrix[,i] <- c(sd_dx[,i - 1], sd_xx[,i - 1])
  }
  
  block_2_matrix <- matrix(c(sd_yd, sd_yx[1,]), ncol = 1)
  
  # Solve for the corrected parameters
  solve_param <- solve(block_1_matrix) %*% block_2_matrix
  solve_param
  
  # Compute the intercept
  intercept <- mean(outcome) -
    ((solve_param[1,1]) * (mean(m_matrix[,1]) - median(pistar12)) / (1 - median(pistar12) - median(pistar21))) -
    t(colMeans(predictor_matrix) %*% solve_param[-1,])
  
  # Intercept not estimated in "solve_param".
  ## Fix this to allow for flexible parameter dimensions! 
  OLS_results <- data.frame(Parameter = c(COMBO_EM_results$Parameter,
                                          "theta0", "theta_m", "theta_x",
                                          paste0("theta_c",
                                                 1:ncol(matrix(c(c_matrix),
                                                               nrow = sample_size,
                                                               byrow = FALSE)))),
                            Estimates = c(c(predicted_beta), c(predicted_gamma),
                                          intercept, solve_param[, 1]),
                            Convergence = rep(COMBO_EM_results$Convergence[1],
                                              n_param),
                            Method = "OLS")
  
  return(OLS_results)
}
