# code_estPS_IR2015
# Estimate PS using logistic regression with variable selection based on 
#   the procedure described by Imbens and Rubin (2015) in Section 13.3.


# Issues ------------------------------------------------------------------
# fn_estPS_IR2015
# [x] Shall we distinguish binary covariates from the other categoricals?
#     - It probably does not matter.
# [x] Shall we include the estimation of PS?
#     - No for now.  Only the terms selected
# [ ] Shall we include a version with the weakly informative prior glm?
# [x] Will there be issues in calculating log-likelihoods (e.g., NA/Inf)?
#     - It seems that the problem was no checking of the empty term_candi.

para_PSmodel_L = 1 # for test statistic in selecting main effects
para_PSmodel_Q = 2.71 # for test statistic in selecting 2nd-order terms

fn_estPS_IR2015 = function(dat, C_L = para_PSmodel_L, C_Q = para_PSmodel_Q) {
  # Select main effects and their 2nd-order terms (squared terms and 2-way)
  #   interaction terms for logistic regression and PS estimation.
  # 
  # Input: 
  #   - data: Data frame with the following columns
  #       - W: Treatment assignment
  #       - Covariates of proper types (continuous or factor)
  # Output: 
  #   - Selected effects
  
  fn_lnLik_new = function(X_new) {
    # Compute the log-likelihood with an additional covariate to the 
    #   current formula.
    # Input : X_new, a string of variable name
    # Output: The log-likelihood
    # Note  : Only used within fn_estPS_IR2015 for code simplicity
    formu_now_elt = as.character(formu_now)
    formu_new = as.formula(paste(formu_now_elt[2], ' ~ ', 
                                 formu_now_elt[3], ' + ', X_new,
                                 sep = ''))
    lnLik_new = logLik(glm(formu_new, data = dat, family = binomial))
    return(lnLik_new)
  }
  
  # The L pass
  
  term_selected = '1'
  term_candi = setdiff(colnames(dat), 'W')
  
  while (TRUE) {
    
    formu_now = as.formula(paste('W ~ ', 
                                 paste(term_selected, collapse = ' + '), 
                                 sep = ''))
    model_now = glm(formu_now, data = dat, family = binomial)
    lnLik_now = logLik(model_now)
    
    lnLik_candi = sapply(term_candi, fn_lnLik_new)
    
    LRT_stat = 2 * (lnLik_candi - lnLik_now)
    
    if (LRT_stat[which.max(LRT_stat)] > C_L) {
      term_new = names(LRT_stat)[which.max(LRT_stat)]
      # print(paste('L pass: Select ', term_new, '.', sep = '')) # for debugging
    } else {
      term_new = NULL
      # print('L pass: No term selected!') # for debugging
      break  
    }
    
    term_selected = union(term_selected, term_new)
    term_candi = setdiff(term_candi, term_new)
    
    if (length(term_candi) == 0) {break}
    
  }
  
  # print(paste('L pass selected: ', paste(term_selected, collapse = ', '), sep = '')) # for debugging
  
  # The Q pass
  
  term_all_fac_ind = sapply(dat, is.factor)
  term_selected_fac = intersect(setdiff(term_selected, '1'), 
                                names(term_all_fac_ind)[term_all_fac_ind == TRUE])
  term_selected_num = setdiff(setdiff(term_selected, '1'),
                              term_selected_fac)
  term_candi_interact = outer(setdiff(term_selected, '1'),
                              setdiff(term_selected, '1'), 
                              paste, sep = ':')
  term_candi = c(paste('I(', term_selected_num, '^2)', sep = ''), 
                 as.vector(term_candi_interact[upper.tri(term_candi_interact)]))
  
  while (TRUE) {
    formu_now = as.formula(paste('W ~ ', 
                                 paste(term_selected, collapse = ' + '), 
                                 sep = ''))
    model_now = glm(formu_now, data = dat, family = binomial)
    lnLik_now = logLik(model_now)

    lnLik_candi = sapply(term_candi, fn_lnLik_new)

    LRT_stat = 2 * (lnLik_candi - lnLik_now)

    if (LRT_stat[which.max(LRT_stat)] > C_Q) {
      term_new = names(LRT_stat[which.max(LRT_stat)])
      # print(paste('Q pass: Select ', term_new, '.', sep = '')) # for debugging
    } else {
      term_new = NULL
      # print('Q pass: No term selected!') # for debugging
      break
    }

    term_selected = union(term_selected, term_new)
    term_candi = setdiff(term_candi, term_new)
    
    if (length(term_candi) == 0) {break}
  
  }
  
  return(term_selected)
  
}