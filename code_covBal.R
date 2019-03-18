# code_covBal
# Some functions for covariate balance checks

# Functions for calculating balance statistics ----------------------------

fn_stdDiff = function(W, cov) {
  # Calculate standardized difference in means
  # Input :
  #   W = binary vector for treatment assignment (TRUE for treated)
  #   cov = vector of covariate values_c
  # Output: numeric, standardized difference in means 
  # Reference: Imben and Rubin (2015)

  ### Under construction
  # [ ] What to do when (xbar_t != xbar_c) & ((ssq_t + ssq_c) == 0)

  
  if ((length(unique(W)) == 1) | 
      (sum(table(W) == 1) > 0)) { # Check for empty/singleton group
    # This only happens to the discarded stratum, thanks to the PSIR algorithm.
    myvalue = 0
  } else {  

    xbar_t = mean(cov[W == TRUE])
    xbar_c = mean(cov[W == FALSE])
    ssq_t = var(cov[W == TRUE])
    ssq_c = var(cov[W == FALSE])

    if ((ssq_t + ssq_c) > 0) {
      myvalue = (xbar_t - xbar_c) / sqrt((ssq_t + ssq_c)/2)
    } else if ((ssq_t + ssq_c) == 0) {
      if (xbar_t == xbar_c) {
        myvalue = 0
      } else if (xbar_t != xbar_c) {
        myvalue = xbar_t - xbar_c
      }
    }
  }
  return(myvalue)
  
}

fn_lnSRatio = function(W, cov) {
  # 
  # Input :
  #   W = binary vector for treatment assignment (TRUE for treated)
  #   cov = vector of covariate values_c
  # Output:
  # Reference: 
  
  ssq_t = var(cov[W == TRUE])
  ssq_c = var(cov[W == FALSE])
  myvalue = log(sqrt(ssq_t / ssq_c))
  return(myvalue)

}

fn_overlap = function(W, cov, alpha = 0.05) {
  #  
  # Input :
  #   W = binary vector for treatment assignment (TRUE for treated)
  #   cov = vector of covariate values_c
  # Output:
  # Reference: 
  
  values_t = cov[W == TRUE]
  values_c = cov[W == FALSE]
  ecdf_t = ecdf(values_t)
  ecdf_c = ecdf(values_c)
  pi_t = 1 - ecdf_t(quantile(values_c, 1 - alpha/2, type = 1)) + ecdf_t(quantile(values_c, alpha/2, type = 1))
  pi_c = 1 - ecdf_c(quantile(values_t, 1 - alpha/2, type = 1)) + ecdf_c(quantile(values_t, alpha/2, type = 1))
  return(c(pi_c = pi_c, pi_t = pi_t))

}