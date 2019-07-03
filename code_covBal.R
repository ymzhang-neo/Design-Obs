# code_covBal
# Some functions for covariate balance checks


# Under construction ------------------------------------------------------
# Tasks to be done are listed here.
# [ ] Review the details and polish fn_lnSRatio, fn_overlap


# Functions for calculating balance statistics ----------------------------

fn_stdDiff = function(W, X_vec) {
  # Calculate standardized difference in means
  # Input :
  #   W = binary vector for treatment assignment (1 for treated)
  #   X_vec = vector of covariate values_c
  # Output: numeric, standardized difference in means 
  # Reference: Imben and Rubin (2015)

  ### Under construction
  # [ ] What to do when (xbar_t != xbar_c) & ((ssq_t + ssq_c) == 0)

  W = as.numeric(W) # in case W is logical
  
  if ((length(unique(W)) == 1) | 
      (sum(table(W) == 1) > 0)) { # Check for empty/singleton group
    # This only happens to the discarded stratum, thanks to the PSIR algorithm.
    myvalue = 0
  } else {  

    xbar_t = mean(X_vec[W == 1], na.rm = TRUE)
    xbar_c = mean(X_vec[W == 0], na.rm = TRUE)
    ssq_t = var(X_vec[W == 1], na.rm = TRUE)
    ssq_c = var(X_vec[W == 0], na.rm = TRUE)

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


fn_lnSRatio = function(W, X_vec) {
  # 
  # Input :
  #   W = binary vector for treatment assignment (TRUE for treated)
  #   X_vec = vector of covariate values_c
  # Output:
  # Reference: 
  
  ssq_t = var(X_vec[W == TRUE])
  ssq_c = var(X_vec[W == FALSE])
  myvalue = log(sqrt(ssq_t / ssq_c))
  return(myvalue)

}


fn_overlap = function(W, X_vec, alpha = 0.05) {
  #  
  # Input :
  #   W = binary vector for treatment assignment (TRUE for treated)
  #   X_vec = vector of covariate values_c
  # Output:
  # Reference: 
  
  values_t = X_vec[W == TRUE]
  values_c = X_vec[W == FALSE]
  ecdf_t = ecdf(values_t)
  ecdf_c = ecdf(values_c)
  pi_t = 1 - ecdf_t(quantile(values_c, 1 - alpha/2, type = 1)) + ecdf_t(quantile(values_c, alpha/2, type = 1))
  pi_c = 1 - ecdf_c(quantile(values_t, 1 - alpha/2, type = 1)) + ecdf_c(quantile(values_t, alpha/2, type = 1))
  return(c(pi_c = pi_c, pi_t = pi_t))

}


fn_stdDiff_subAgg = function(W, X_vec, stra, by_trt = TRUE) {
  # Compute the aggregated stdDiff of given stratification
  # Input :
  #   W = binary vector for treatment assignment (TRUE for treated)
  #   X_vec = vector of covariate values_c
  #   stra = vector of integer-valued stratum labels
  #   by_trt = logical, aggregation by treated group sizes if TRUE,
  #     aggregation by stratum sizes if FALSE
  # Output: numeric, stdDiff in means aggregated over strata
  
  # Under construction
  # [ ] What to do when stra == 0? [temp] Ignore stra == 0
  
  W = W[stra != 0]
  X_vec = X_vec[stra != 0]
  stra = stra[stra != 0]
  
  dat_split = split(data.frame(W=W, X_vec=X_vec), stra)
  
  stdDiff_byStra = sapply(dat_split, function(oneStra) {
    fn_stdDiff(W = oneStra$W, X_vec = oneStra$X_vec)
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



# Visualization -----------------------------------------------------------


fn_covbal_hist = function(W, X_vec, 
                          probability=TRUE, overlay=FALSE, 
                          xlab = 'Covariate', 
                          main = 'Comparative histograms',
                          include_legend = TRUE,
                          col = para_col) {
  # Comparative histograms for a continuous covariate.
  
  plot_W0 = hist(X_vec[W == 0], plot = FALSE)
  plot_W1 = hist(X_vec[W == 1], plot = FALSE)
  
  plot_xlim = range(X_vec)
  if (probability) {
    plot_ylim = c(0, max(max(plot_W0$density), max(plot_W1$density)))
  } else {
    plot_ylim = c(0, max(max(plot_W0$counts), max(plot_W1$counts)))
  }
  
  if (overlay) {
    hist(X_vec[W == 0], probability = probability, col = col[1], 
         xlim = plot_xlim, ylim = plot_ylim, xlab = xlab, main = main)
    hist(X_vec[W == 1], probability = probability, col = col[2], 
         xlim = plot_xlim, ylim = plot_ylim, xlab = xlab, add = TRUE)
  } else {
    par(mfrow = c(1, 2))
    hist(X_vec[W == 0], probability = probability, col = col[1], 
         xlim = plot_xlim, ylim = plot_ylim, xlab = xlab, main = 'W == 0')
    hist(X_vec[W == 1], probability = probability, col = col[2], 
         xlim = plot_xlim, ylim = plot_ylim, xlab = xlab, main = 'W == 1')
  }
  
  if (include_legend) {
    legend('topright', c('Control', 'Treated'), pch = 15, col = col)
  }
  
  par(mfrow = c(1, 1))
  
}
