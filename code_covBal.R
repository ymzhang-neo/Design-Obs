# code_covBal
# Some functions for covariate balance checks


# Under construction ------------------------------------------------------
# Tasks to be done are listed here.
# [ ] Review the details and polish fn_lnSRatio, fn_overlap


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

fn_stdDiff_subAgg = function(W, cov, stra, by_trt = TRUE) {
  # Compute the aggregated stdDiff of given stratification
  # Input :
  #   W = binary vector for treatment assignment (TRUE for treated)
  #   cov = vector of covariate values_c
  #   stra = vector of integer-valued stratum labels
  #   by_trt = logical, aggregation by treated group sizes if TRUE,
  #     aggregation by stratum sizes if FALSE
  # Output: numeric, stdDiff in means aggregated over strata
  
  # Under construction
  # [ ] What to do when stra == 0? [temp] Ignore stra == 0
  
  W = W[stra != 0]
  cov = cov[stra != 0]
  stra = stra[stra != 0]
  
  dat_split = split(data.frame(W=W, cov=cov), stra)
  
  stdDiff_byStra = sapply(dat_split, function(oneStra) {
    fn_stdDiff(W = oneStra$W, cov = oneStra$cov)
  })

  if (by_trt) {
    straWt_num = sapply(dat_split, function(oneStra) {sum(oneStra$W)})
    straWt = straWt_num / sum(straWt_num)
  } else {
    straWt_num = sapply(dat_split, function(oneStra) {length(oneStra$W)})
    straWt = straWt_num / sum(straWt_num)
  }
  
  stdDiff_agg = sum(stdDiff_byStra * straWt)
  
  return(list(stdDiff_byStra = stdDiff_byStra, 
              straWt = straWt, 
              stdDiff_agg = stdDiff_agg))
  
}
