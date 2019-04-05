# code_estTrt_Neyman
# Estimate treatment effect using the Neyman approach

fn_estTrt_Neyman = function(Y, W, stra = NULL, clevel = 0.95) {
  # Input :
  #   Y = the vector of response value, assumed continuous for now
  #   W = the vector of treatment assignment (1 for treated)
  #   stra = the vector of subgroup labels (0 for discarded)
  # Output: point estimate and confidence interval
  # 
  # Under construction
  # - Incorporate handling of binary Y  
  # - Different estimator of sampling variance
  # - Weighting choice: treated or stratum sizes
  
  if ((is.null(stra)) | (length(setdiff(unique(stra), 0)) == 1)) {
    if (is.null(stra)) {
      estTrt = mean(Y[W == 1], na.rm = TRUE) - mean(Y[W == 0], na.rm = TRUE)
      estTrt_ssq = var(Y[W == 1], na.rm = TRUE) / (sum(W == 1) - 1) + 
        var(Y[W == 0], na.rm = TRUE) / (sum(W == 0) - 1)  
    } else if (length(setdiff(unique(stra), 0)) == 1) {
      Y = Y[stra != 0]
      W = W[stra != 0]
      estTrt = mean(Y[W == 1], na.rm = TRUE) - mean(Y[W == 0], na.rm = TRUE)
      estTrt_ssq = var(Y[W == 1], na.rm = TRUE) / (sum(W == 1) - 1) + 
        var(Y[W == 0], na.rm = TRUE) / (sum(W == 0) - 1)  
    }
  } else {
    dat_split = split(data.frame(Y = Y, W = W), stra)
    estTrt_base = sapply(dat_split, function(oneStra) {
      if (sum(oneStra$W == 1) * sum(oneStra$W == 0) != 0) {
        oneEst = mean(oneStra$Y[oneStra$W == 1], na.rm = TRUE) - 
          mean(oneStra$Y[oneStra$W == 0], na.rm = TRUE)
        oneEst_ssq = var(oneStra$Y[oneStra$W == 1], na.rm = TRUE) / (sum(oneStra$W == 1) - 1) + 
          var(oneStra$Y[oneStra$W == 0], na.rm = TRUE) / (sum(oneStra$W == 0) - 1)
      } else {
        oneEst = 0
        oneEst_ssq = 0
      }
      return(c(est=oneEst, ssq=oneEst_ssq))
    })
    
    tab_W_stra = table(W, stra)
    if (0 %in% stra) {
      stra_wt = c(`0` = 0, tab_W_stra[2, ][-1] / sum(tab_W_stra[2, ][-1]))
    } else {
      stra_wt = tab_W_stra[2, ] / sum(tab_W_stra[2, ])
    }
    
    estTrt = sum(estTrt_base[1, ] * stra_wt)
    estTrt_ssq = sum(estTrt_base[2, ] * ((stra_wt)^2))
    
  }
  
  CI_coef = qnorm((1 - clevel)/2, lower.tail = FALSE)
  CI_low = estTrt - CI_coef * sqrt(estTrt_ssq)
  CI_high = estTrt + CI_coef * sqrt(estTrt_ssq)
  
  return(list(estTrt = estTrt, CI = c(CI_low, CI_high)))
  
}
