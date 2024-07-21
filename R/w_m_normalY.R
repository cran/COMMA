#' Compute E-step for Binary Mediator Misclassification Model Estimated With the EM Algorithm
#' 
#' Note that this function should only be used for Normal outcome models.
#'
#' @param mstar_matrix A numeric matrix of indicator variables (0, 1) for the observed
#'   mediator \code{M*}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed mediator category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param pistar_matrix A numeric matrix of conditional probabilities obtained from
#'   the internal function \code{pistar_compute}. Rows of the matrix correspond to
#'   each subject and to each observed mediator category. Columns of the matrix
#'   correspond to each true, latent mediator category.
#' @param pi_matrix A numeric matrix of probabilities obtained from the internal
#'   function \code{pi_compute}. Rows of the matrix correspond to each subject.
#'   Columns of the matrix correspond to each true, latent mediator category.
#' @param p_yi_m0 A numeric vector of Normal outcome likelihoods computed assuming a
#'   true mediator value of 0.
#' @param p_yi_m1 A numeric vector of Normal outcome likelihoods computed assuming a
#'   true mediator value of 1.
#' @param sample_size An integer value specifying the number of observations in
#'   the sample. This value should be equal to the number of rows of the observed
#'   mediator matrix, \code{mstar_matrix}.
#' @param n_cat The number of categorical values that the true outcome, \code{M},
#'   and the observed outcome, \code{M*}, can take.
#'
#' @return \code{w_m_normalY} returns a matrix of E-step weights for the EM-algorithm.
#'   Rows of the matrix correspond to each subject. Columns of the matrix correspond
#'   to the true mediator categories \eqn{j = 1, \dots,} \code{n_cat}.
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#'
#' @importFrom stats rnorm rgamma rmultinom
w_m_normalY <- function(mstar_matrix, # Observed mediator matrix
                        pistar_matrix, # Probability of observed mediator given latent mediator
                        pi_matrix, # Mediator probabilities
                        p_yi_m0, p_yi_m1, # Likelihood value of Y for given M
                        sample_size, n_cat){
  
  y_m1_mstar1_frac1 = p_yi_m1 * mstar_matrix[,1] * pistar_matrix[1:sample_size, 1] *
    pi_matrix[,1]
  y_m2_mstar1_frac1 = p_yi_m0 * mstar_matrix[,1] * pistar_matrix[1:sample_size, 2] *
    pi_matrix[,2]
  
  y_m1_mstar2_frac1 = p_yi_m1 * mstar_matrix[,2] *
    pistar_matrix[(sample_size + 1):(n_cat * sample_size), 1] * pi_matrix[,1]
  y_m2_mstar2_frac1 = p_yi_m0 * mstar_matrix[,2] *
    pistar_matrix[(sample_size + 1):(n_cat * sample_size), 2] * pi_matrix[,2]
  
  mstar1_frac2 = (p_yi_m1 * pistar_matrix[1:sample_size, 1] * pi_matrix[,1]) +
    (p_yi_m0 * pistar_matrix[1:sample_size, 2] * pi_matrix[,2])
  
  mstar2_frac2 = (p_yi_m0 *
                    pistar_matrix[(sample_size + 1):(n_cat * sample_size), 2]
                  * pi_matrix[,2]) +
    (p_yi_m1 * pistar_matrix[(sample_size + 1):(n_cat * sample_size), 1] * pi_matrix[,1])
  
  m1_term = (y_m1_mstar1_frac1 / mstar1_frac2) + (y_m1_mstar2_frac1 / mstar2_frac2)
  m2_term = (y_m2_mstar1_frac1 / mstar1_frac2) + (y_m2_mstar2_frac1 / mstar2_frac2)
  
  weight_m = matrix(c(m1_term, m2_term), nrow = sample_size, byrow = FALSE)
  
  return(weight_m)
}
