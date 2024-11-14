## Group members : 1.Yitian Zhu 2.Duoer Chen 3.Yaqi Liu
## Contributions : 1.Yitian Zhu 34% 2.Duoer Chen 34% 3.Yaqi Liu 32%
## Explanations : We roughly divided the work as follows:
## Yitian Zhu was mainly responsible for writing the function LMMprof.
## Duoer Chen was responsible for writing the function LMMsetup.
## Yaqi Liu was responsible for writing the function lmm and modifying the code.
## But throughout the work we worked together on the project.
## After extensive discussion, debugging and polishing,
## we finally arrived at the current outcome.



# Main function: lmm
# This function performs linear mixed model (LMM) estimation.
# It estimates the fixed effects (beta) and random effects (theta) using maximum likelihood estimation.
# Arguments:
# - form: Model formula specifying fixed effects.
# - dat: Data frame containing the data to be analyzed.
# - ref: List of variables specifying random effects (optional).
lmm <- function(form, dat, ref = list()) {
  # Step 1: Set up model matrices X  and Z, and initial parameter values for theta
  setup <- LMMsetup(form, dat, ref)
  
  # store the final estimate of beta_hat 
  final_beta_hat <- NULL
  
  # Step 2: Use the optim function to maximize the likelihood with respect to theta 
  opt_result <- optim(
    par = setup$theta_init,  # Initial theta values 
    fn = function(par) {
      # Calculate the log-likelihood for the current theta values
      result <- LMMprof(par, setup = setup)
      
      # Extract and store beta_hat as an external variable
      final_beta_hat <<- attr(result, "beta_hat")  
      
      return(result)  # Return the computed log-likelihood
    }
  )
  
  # Step 3: Return the estimated parameters
  list(beta = final_beta_hat, theta = opt_result$par)
}




# Helper function: LMMsetup
# This function sets up the design matrices (X, Z) and initializes theta.
# It returns a list containing all the necessary setup for the LMMprof estimation.
# Arguments:
# - form: Model formula (specifies the fixed effects).
# - dat: Data frame containing the data.
# - ref: List of variables for random effects (optional).
LMMsetup <- function(form, dat, ref) {
  # Extract the name of the response variable from the formula 
  response_var <- as.character(form[[2]])
  
  # Convert the response variable data into a numeric vector (y)
  y <- as.numeric(dat[[response_var]])
  
  # Create the fixed effects design matrix X based on the formula
  X <- model.matrix(form, data = dat)
  
  if (length(ref) == 0) {
    # Case with no random effects:
    Z <- NULL
    theta_init <- 1  # Initialize with only residual variance
    random_effect_col_counts <- NULL
  } 
  
  else {
    # Case with random effects:
    # Generate Z matrices for each random effect variable
    Z_matrices <- lapply(ref, function(vars) {
      model.matrix(as.formula(paste("~", paste(vars, collapse = ":"), "- 1")), data = dat)
    })
    
    # Combine all Z matrices by columns to form the full Z matrix
    Z <- do.call(cbind, Z_matrices)
    
    # Record column count for each random effect to handle multiple random effects
    random_effect_col_counts <- sapply(Z_matrices, ncol)
    
    # Initialize theta with log-transformed random effect variances and residual variance
    theta_init <- rep(log(0.5), length(ref) + 1)
  }
  
  
  list(X = X, Z = Z, y = y, theta_init = theta_init, col = random_effect_col_counts)
}



# Profile log-likelihood function: LMMprof
# This function calculates the negative log-likelihood given the current values of theta.
# It uses the setup information to compute the log-likelihood for a linear mixed model.
# Arguments:
# - theta: The current parameter values (residual variance and random effects variances).
# - setup: A list containing the setup information (design matrices and needed information).
LMMprof <- function(theta, setup) {
  X <- setup$X
  Z <- setup$Z
  y <- setup$y
  n <- length(y)  
  col <- setup$col
  p <- ncol(Z)
  
  # Decompose theta into residual variance (sigma2) and random effect variances
  sigma2 <- exp(2 * theta[1])  
  
  # Case with no random effects:
  if (is.null(Z)) {
    # Perform ordinary least squares (OLS) estimation for beta (fixed effects)
    cholx <- chol(t(X) %*% X)
    xty <- t(X) %*% y
    beta_hat <- backsolve(cholx, forwardsolve(t(cholx), xty))
    
    # Calculate residuals
    residual <- y - X %*% beta_hat
    rWr <- sum(residual^2) / sigma2  # Residual sum of squares divided by variance
    
    # Total negative log-likelihood for OLS case
    term1 <- -0.5 * rWr  
    term2 <- -0.5 * n * log(sigma2)
    log_likelihood <- -(term1 + term2)  
  } 
  
  
  # Case with random effects:
  else {
    random_vars <- exp(2 * theta[-1])  # Random effect variances
    psi_theta <- diag(rep(random_vars, col), p, p)
    
    # QR decomposition of Z to obtain R matrix for stability
    qr_Z <- qr(Z)
    R <- qr.R(qr_Z)
    
    # Construct A = RψθR^T + Iσ^2 
    A <- R %*% psi_theta %*% t(R) + diag(sigma2, p, p)
    chol_A <- chol(A)  
    
    # Compute (Q^T * y and Q^T * X) for further calculations
    Qty <- qr.qty(qr_Z, y)
    Qtx <- qr.qty(qr_Z, X)
    
    # Compute W's diagonal terms (using Cholesky decomposition)
    upper_block <- backsolve(chol_A, forwardsolve(t(chol_A), diag(1, p, p)))
    lower_block <- diag(1 / sigma2, n - p, n - p)
    diag_W <- as.matrix(Matrix::bdiag(upper_block, lower_block))
    
    # Compute Wy and WX efficiently using matrix operations
    Wy <- qr.qy(qr_Z, diag_W %*% Qty)
    Wx <- qr.qy(qr_Z, diag_W %*% Qtx)
    
    # Compute the necessary cross-products x^T W y and x^T W x
    xWy <- t(X) %*% Wy
    xWx <- t(X) %*% Wx
    
    # Cholesky decomposition to solve for beta_hat
    chol_xWx <- chol(xWx)
    beta_hat <- backsolve(chol_xWx, forwardsolve(t(chol_xWx), xWy))
    
    # Calculate residuals and additional log-likelihood components
    residual <- y - X %*% beta_hat
    Qtr <- qr.qty(qr_Z, residual)
    Wr <- qr.qy(qr_Z, diag_W %*% Qtr)
    rWr <- t(residual) %*% Wr
    
    # Compute the total negative log-likelihood
    term1 <- -0.5 * rWr
    log_det_A <- sum(log(diag(chol_A))) * 2  
    log_det_sigma <- (n - p) * log(sigma2)  
    term2 <- -0.5 * (log_det_A + log_det_sigma)  
    
    log_likelihood <- -(term1 + term2)
  }
  
  # Store the estimated beta_hat as an attribute of the log-likelihood
  attr(log_likelihood, "beta_hat") <- beta_hat
  
  return(log_likelihood)  
}




