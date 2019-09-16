# code_covBal
# Some functions for covariate balance checks


# Under construction ------------------------------------------------------
# Tasks to be done are listed here.
# [ ] Review the details and polish fn_lnSRatio, fn_overlap


# Functions for calculating balance statistics ----------------------------

fn_stdDiff = function(W, X_vec, use_abs = FALSE) {
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
  
  if (use_abs) {
    myvalue = abs(myvalue)
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

fn_stdDiff_overall = function(covbal) {
  # Compute an overall balance measure based on those of all covariates.
  
  # covbal_overall = sum(abs(covbal)) / length(covbal)
  covbal_overall = sqrt(sum(covbal^2) / length(covbal)) # larger weights for covariates with more severe imbalance
  return(covbal_overall)
}


# Visualization -----------------------------------------------------------


fn_covbal_hist = function(W, X_vec, 
                          probability = TRUE, overlay = FALSE, 
                          xlab = 'Covariate', 
                          main = 'Comparative histograms',
                          include_legend = TRUE,
                          col = para_col) {
  # Comparative histograms for a continuous covariate.
  # [TODO] Allow outside change of plot panel 
  #   Currently overlay= would determine "par(mfrow=))"
  
  para_col = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2))
  
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
    par_current = par(no.readonly = TRUE)
    par(mfrow = c(1, 2))
    hist(X_vec[W == 0], probability = probability, col = col[1], 
         xlim = plot_xlim, ylim = plot_ylim, xlab = xlab, main = 'W == 0')
    hist(X_vec[W == 1], probability = probability, col = col[2], 
         xlim = plot_xlim, ylim = plot_ylim, xlab = xlab, main = 'W == 1')
    par(par_current)
  }
  
  if (include_legend) {
    legend('topright', c('Control', 'Treated'), pch = 15, col = col)
  }
  
}

fn_covbal_scatter2d = function(W, X_vec, Y_vec, stra, 
                               label_stra = TRUE, col = para_col) {
  # Scatter plot of two covariates, colored by the treatment status, 
  #   labeled by the strata.
  
  if (label_stra) {
    plot(X_vec, Y_vec, col = 'white', pch = 16)
    text(X_vec, Y_vec, labels = stra, col = para_col[W + 1])
  } else {
    plot(X_vec, Y_vec, col = para_col[W + 1], pch = 16)
  }
  
}

fn_vis_stdDiff_comp0 = function(covbal1, covbal2, 
                               xlim=NULL, xlab=NULL, main=NULL, 
                               horizontal_ref=FALSE,
                               vertical_ref=c(-0.1, 0.1),
                               arrow=TRUE, arrow_length=0.1, 
                               name1=NULL, name2=NULL, 
                               decimals=3) {
  # Compare standardized differences by two subclassifications
  # For empty imput of covbal2, use covbal2=NA
  
  n_covariate = length(covbal1)
  
  # Temporary: covbal2=NA means no plot for the second design.
  if (length(covbal2) == 1) {
    if (is.na(covbal2)) {
      covbal2 = rep(NA, n_covariate)
    }
  }
  
  # result_col = c('orangered', 'mediumseagreen')
  # result_col_id = (abs(covbal2) <= abs(covbal1)) + 1
  result_col_list = c('orangered', 'mediumseagreen', 
                      'lightcoral')
  result_col = rep('mediumseagreen', length(covbal1))
  result_col[(abs(covbal2) > abs(covbal1)) & 
               (abs(covbal2) <= vertical_ref[2])] = 'lightcoral'
  result_col[(abs(covbal2) > abs(covbal1)) & 
               (abs(covbal2) > vertical_ref[2])] = 'orangered'
  
  if (is.null(xlim)) {
    xlim = max(range(covbal1, covbal2, na.rm = TRUE)) * c(-1, 1)
  }
  
  if (is.null(xlab)) {
    xlab = 'Standardized difference'
  }
  
  if (is.null(main)) {
    main = 'Change of covariate balance'  
  }
  
  # yticks_at = seq(n_covariate, 1, length.out = 10)
  yticks_at = seq(n_covariate, 1, by = -10)
  if (!is.null(names(covbal1))) {
    # yticks = names(covbal1)[seq(1, n_covariate, length.out = 10)]
	yticks = names(covbal1)[seq(1, n_covariate, by = 10)]
  } else {
    # yticks = paste0('X_', seq(1, n_covariate, length.out = 10))
	yticks = paste0('X_', seq(1, n_covariate, by = 10))
  }
  
  plot(x = covbal1, y = n_covariate:1, pch = 16, col = 'gray', 
       xlim = xlim, 
       xlab = xlab,
       ylab = 'Covariate', yaxt = 'n', 
       main = main)
  points(x = covbal2, y = n_covariate:1, pch = 16, 
         col = result_col) # result_col[result_col_id])
  
  if (arrow) {
    arrows(x0 = covbal1, y0 = n_covariate:1, x1 = covbal2, 
           length = arrow_length, col = 'gray35')
  } else {
    segments(x0 = covbal1, y0 = n_covariate:1, x1 = covbal2, 
	         col = 'gray35')
  }
  
  axis(2, at = yticks_at, labels = yticks, las = 1)
  abline(v = 0, col = 'gray')
  abline(v = vertical_ref, col = 'gray', lty = 2)
  
  if (horizontal_ref) {
    abline(h = n_covariate:1, col = 'gray', lty = 2)
  }
  
  legend('bottomright', c('Before', 'Better', 'Worse'), 
         pch = 16, col = c('gray', 'mediumseagreen', 'orangered'))
  
  covbal1_overall = fn_stdDiff_overall(covbal1)
  covbal2_overall = fn_stdDiff_overall(covbal2)
  title(sub = paste0('Overall balance: ', 
                     'Before ', round(covbal1_overall, decimals), ', ', 
                     'After ', round(covbal2_overall, decimals)),
        adj = 1, line = 2, font = 1)
  
}


fn_vis_stdDiff_comp = function(dat, W=NULL, stra1, stra2, 
                               xlim=NULL, xlab=NULL, main=NULL, 
                               horizontal_ref=FALSE,
                               vertical_ref=c(-0.1, 0.1),
                               arrow=TRUE, arrow_length=0.1, 
                               name1=NULL, name2=NULL,
                               decimals=3, add_legend=FALSE) {
  # Compare standardized differences by two subclassifications
  # 
  # Parameters
  # ----------
  # W : NULL or character, optional (default=NULL)
  #   The column name of treatment status in dat.
  #   - If W=NULL, it is assumed that dat has a column 'W' for treatment.
  #   - If W is character, it will be used as the name of a column in dat.
  # 
  # Notes
  # -----
  # - This is an alternative version that allow users to use the raw data.
  # [TODO] Allow one stratification proposal
  
  if (is.null(W)) {
    dat_X = dat[, setdiff(colnames(dat), 'W')]
	W = dat[, 'W']
  } else if (is.character(W)) {
    dat_X = dat[, setdiff(colnames(dat), W)]
	W = dat[, W]
  } 
  
  n_covariate = ncol(dat_X)
  
  covbal1 = apply(dat_X, 2, function(one_column){
    return(fn_stdDiff_subAgg(W = W, X_vec = one_column, stra = stra1)$stdDiff_agg)
  })
  covbal2 = apply(dat_X, 2, function(one_column){
    return(fn_stdDiff_subAgg(W = W, X_vec = one_column, stra = stra2)$stdDiff_agg)
  })
  
  covbal1_overall = fn_stdDiff_overall(covbal1)
  covbal2_overall = fn_stdDiff_overall(covbal2)
  
  keep_n1 = sum(stra1 != 0)
  keep_n2 = sum(stra2 != 0)
  discard_n1 = sum(stra1 == 0)
  discard_n2 = sum(stra2 == 0)
  total_n1 = length(stra1) 
  total_n2 = length(stra2)
  
  #[TODO] Summary over all covariates
  
  result_col_list = c('orangered', 'mediumseagreen', 
                 'lightcoral')
  result_col = rep('mediumseagreen', length(covbal1))
  result_col[(abs(covbal2) > abs(covbal1)) & 
               (abs(covbal2) <= vertical_ref[2])] = 'lightcoral'
  result_col[(abs(covbal2) > abs(covbal1)) & 
               (abs(covbal2) > vertical_ref[2])] = 'orangered'
  
  if (is.null(xlim)) {
    xlim = max(range(covbal1, covbal2, na.rm = TRUE)) * c(-1, 1)
  }
  
  if (is.null(xlab)) {
    xlab = 'Standardized difference'
  }
  
  if (is.null(main)) {
    main = 'Change of covariate balance'  
  }
  
  # yticks_at = seq(n_covariate, 1, length.out = 10)
  yticks_at = seq(n_covariate, 1, by = -10)
  if (!is.null(names(covbal1))) {
    # yticks = names(covbal1)[seq(1, n_covariate, length.out = 10)]
	yticks = names(covbal1)[seq(1, n_covariate, by = 10)]
  } else {
    # yticks = paste0('X_', seq(1, n_covariate, length.out = 10))
	yticks = paste0('X_', seq(1, n_covariate, by = 10))
  }
  
  plot(x = covbal1, y = n_covariate:1, pch = 16, col = 'gray', 
       xlim = xlim, 
       xlab = xlab,
       ylab = 'Covariate', yaxt = 'n', 
       main = main)
  points(x = covbal2, y = n_covariate:1, pch = 16, 
         col = result_col) # result_col[result_col_id])
  
  if (arrow) {
    arrows(x0 = covbal1, y0 = n_covariate:1, x1 = covbal2, 
           length = arrow_length, col = 'gray35')
  } else {
    segments(x0 = covbal1, y0 = n_covariate:1, x1 = covbal2, 
	         col = 'gray35')
  }
  
  axis(2, at = yticks_at, labels = yticks, las = 1)
  abline(v = 0, col = 'gray')
  abline(v = vertical_ref, col = 'gray', lty = 2)
  
  if (horizontal_ref) {
    abline(h = n_covariate:1, col = 'gray', lty = 2)
  }
  
  if (add_legend) {
    legend('topright', # 'bottomright', 
           c('Before', 'Better', 'Acceptable', 'Worse'), 
           pch = 16, col = c('gray', 'mediumseagreen', 'lightcoral', 'orangered'))
  }
  
  title(sub = paste0('Sample sizes: ', 
                     'Before ', keep_n1, ', ',
					 'After ', keep_n2), 
        adj = 1, line = 2, font = 1)
  title(sub = paste0('Overall balance: ', 
                     'Before ', round(covbal1_overall, decimals), ', ', 
					 'After ', round(covbal2_overall, decimals)),
        adj = 1, line = 3, font = 1)
  
}


fn_vis_stdDiff = function(dat, W=NULL, stra1, # stra2, 
                          xlim=NULL, xlab=NULL, main=NULL, 
                          horizontal_ref=FALSE,
                          vertical_ref=c(-0.1, 0.1),
                          # arrow=TRUE, arrow_length=0.1, 
                          name1=NULL, # name2=NULL,
                          decimals=3 #,  add_legend=FALSE
) {
  # Standardized differences by two subclassifications, no comparison
  # 
  # Parameters
  # ----------
  # W : NULL or character, optional (default=NULL)
  #   The column name of treatment status in dat.
  #   - If W=NULL, it is assumed that dat has a column 'W' for treatment.
  #   - If W is character, it will be used as the name of a column in dat.
  # 
  # Notes
  # -----
  # - This is an alternative version that allow users to use the raw data.
  # [TODO] Allow one stratification proposal
  
  if (is.null(W)) {
    dat_X = dat[, setdiff(colnames(dat), 'W')]
  W = dat[, 'W']
  } else if (is.character(W)) {
    dat_X = dat[, setdiff(colnames(dat), W)]
  W = dat[, W]
  } 
  
  n_covariate = ncol(dat_X)
  
  covbal1 = apply(dat_X, 2, function(one_column){
    return(fn_stdDiff_subAgg(W = W, X_vec = one_column, stra = stra1)$stdDiff_agg)
  })
  # covbal2 = apply(dat_X, 2, function(one_column){
  #   return(fn_stdDiff_subAgg(W = W, X_vec = one_column, stra = stra2)$stdDiff_agg)
  # })
  
  covbal1_overall = fn_stdDiff_overall(covbal1)
  # covbal2_overall = fn_stdDiff_overall(covbal2)
  
  keep_n1 = sum(stra1 != 0)
  # keep_n2 = sum(stra2 != 0)
  discard_n1 = sum(stra1 == 0)
  # discard_n2 = sum(stra2 == 0)
  total_n1 = length(stra1) 
  # total_n2 = length(stra2)
  
  #[TODO] Summary over all covariates
  
  # result_col_list = c('orangered', 'mediumseagreen', 
  #                'lightcoral')
  # result_col = rep('mediumseagreen', length(covbal1))
  # result_col[(abs(covbal2) > abs(covbal1)) & 
  #              (abs(covbal2) <= vertical_ref[2])] = 'lightcoral'
  # result_col[(abs(covbal2) > abs(covbal1)) & 
  #              (abs(covbal2) > vertical_ref[2])] = 'orangered'
  
  if (is.null(xlim)) {
    xlim = max(range(covbal1, na.rm = TRUE)) * c(-1, 1)
    # xlim = max(range(covbal1, covbal2, na.rm = TRUE)) * c(-1, 1)
  }
  
  if (is.null(xlab)) {
    xlab = 'Standardized difference'
  }
  
  if (is.null(main)) {
    main = 'Covariate balance'  
  }
  
  # yticks_at = seq(n_covariate, 1, length.out = 10)
  yticks_at = seq(n_covariate, 1, by = -10)
  if (!is.null(names(covbal1))) {
    # yticks = names(covbal1)[seq(1, n_covariate, length.out = 10)]
  yticks = names(covbal1)[seq(1, n_covariate, by = 10)]
  } else {
    # yticks = paste0('X_', seq(1, n_covariate, length.out = 10))
  yticks = paste0('X_', seq(1, n_covariate, by = 10))
  }
  
  plot(x = covbal1, y = n_covariate:1, pch = 16, # col = 'gray', 
       xlim = xlim, 
       xlab = xlab,
       ylab = 'Covariate', yaxt = 'n', 
       main = main)
  # points(x = covbal2, y = n_covariate:1, pch = 16, 
  #        col = result_col) # result_col[result_col_id])
  
  # if (arrow) {
  #   arrows(x0 = covbal1, y0 = n_covariate:1, x1 = covbal2, 
  #          length = arrow_length, col = 'gray35')
  # } else {
  #   segments(x0 = covbal1, y0 = n_covariate:1, x1 = covbal2, 
  #          col = 'gray35')
  # }
  
  axis(2, at = yticks_at, labels = yticks, las = 1)
  abline(v = 0, col = 'gray')
  abline(v = vertical_ref, col = 'gray', lty = 2)
  
  if (horizontal_ref) {
    abline(h = n_covariate:1, col = 'gray', lty = 2)
  }
  
  # if (add_legend) {
  #   legend('topright', # 'bottomright', 
  #          c('Before', 'Better', 'Acceptable', 'Worse'), 
  #          pch = 16, col = c('gray', 'mediumseagreen', 'lightcoral', 'orangered'))
  # }
  
  title(sub = paste0('Sample sizes: ', keep_n1), adj = 1, line = 2, font = 1)
  title(sub = paste0('Overall balance: ', round(covbal1_overall, decimals)),
        adj = 1, line = 3, font = 1)
  
}


fn_vis_stdDiff_comp_wide = function(dat, W=NULL, stra1, stra2, 
                                    ylim=NULL, ylab=NULL, main=NULL, 
                                    vertical_ref=FALSE, 
                                    horizontal_ref=c(-0.1, 0.1),
                                    arrow=TRUE, arrow_length=0.1, 
                                    name1=NULL, name2=NULL,
                                    decimals=3, add_legend=FALSE) {
  # Compare standardized differences by two subclassifications
  # 
  # Parameters
  # ----------
  # W : NULL or character, optional (default=NULL)
  #   The column name of treatment status in dat.
  #   - If W=NULL, it is assumed that dat has a column 'W' for treatment.
  #   - If W is character, it will be used as the name of a column in dat.
  # 
  # Notes
  # -----
  # - This is an alternative version that allow users to use the raw data.
  # [TODO] Allow one stratification proposal
  
  if (is.null(W)) {
    dat_X = dat[, setdiff(colnames(dat), 'W')]
  W = dat[, 'W']
  } else if (is.character(W)) {
    dat_X = dat[, setdiff(colnames(dat), W)]
  W = dat[, W]
  } 
  
  n_covariate = ncol(dat_X)
  
  covbal1 = apply(dat_X, 2, function(one_column){
    return(fn_stdDiff_subAgg(W = W, X_vec = one_column, stra = stra1)$stdDiff_agg)
  })
  covbal2 = apply(dat_X, 2, function(one_column){
    return(fn_stdDiff_subAgg(W = W, X_vec = one_column, stra = stra2)$stdDiff_agg)
  })
  
  covbal1_overall = fn_stdDiff_overall(covbal1)
  covbal2_overall = fn_stdDiff_overall(covbal2)
  
  keep_n1 = sum(stra1 != 0)
  keep_n2 = sum(stra2 != 0)
  discard_n1 = sum(stra1 == 0)
  discard_n2 = sum(stra2 == 0)
  total_n1 = length(stra1) 
  total_n2 = length(stra2)
  
  #[TODO] Summary over all covariates
  
  result_col_list = c('orangered', 'mediumseagreen', 
                 'lightcoral')
  result_col = rep('mediumseagreen', length(covbal1))
  result_col[(abs(covbal2) > abs(covbal1)) & 
               (abs(covbal2) <= horizontal_ref[2])] = 'lightcoral'
  result_col[(abs(covbal2) > abs(covbal1)) & 
               (abs(covbal2) > horizontal_ref[2])] = 'orangered'
  
  if (is.null(ylim)) {
    ylim = max(range(covbal1, covbal2, na.rm = TRUE)) * c(-1, 1)
  }
  
  if (is.null(ylab)) {
    ylab = 'Standardized difference'
  }
  
  if (is.null(main)) {
    main = 'Change of covariate balance'  
  }
  
  # yticks_at = seq(n_covariate, 1, length.out = 10)
  xticks_at = seq(1, n_covariate, by = 10)
  if (!is.null(names(covbal1))) {
    # yticks = names(covbal1)[seq(1, n_covariate, length.out = 10)]
    xticks = names(covbal1)[seq(1, n_covariate, by = 10)]
  } else {
    # yticks = paste0('X_', seq(1, n_covariate, length.out = 10))
    xticks = paste0('X_', seq(1, n_covariate, by = 10))
  }
  
  plot(y = covbal1, x = 1:n_covariate, pch = 16, col = 'gray', 
       ylim = ylim, 
       ylab = ylab,
       xaxt = 'n',  xlab = '', #'Covariate', 
       main = main)
  points(y = covbal2, x = 1:n_covariate, pch = 16, 
         col = result_col) # result_col[result_col_id])
  
  if (arrow) {
    arrows(y0 = covbal1, x0 = 1:n_covariate, y1 = covbal2, 
           length = arrow_length, col = 'gray35')
  } else {
    segments(y0 = covbal1, x0 = 1:n_covariate, y1 = covbal2, 
           col = 'gray35')
  }
  
  axis(1, at = xticks_at, labels = xticks, las = 1)
  abline(h = 0, col = 'gray')
  abline(h = horizontal_ref, col = 'gray', lty = 2)
  
  if (vertical_ref) {
    abline(v = n_covariate:1, col = 'gray', lty = 2)
  }
  
  if (add_legend) {
    legend('topright', # 'bottomright', 
           c('Before', 'Better', 'Acceptable', 'Worse'), 
           pch = 16, col = c('gray', 'mediumseagreen', 'lightcoral', 'orangered'))
  }
  
  title(sub = paste0('Sample sizes: ', 
                     'Before ', keep_n1, ', ',
           'After ', keep_n2), 
        adj = 1, line = 2, font = 1)
  title(sub = paste0('Overall balance: ', 
                     'Before ', round(covbal1_overall, decimals), ', ', 
           'After ', round(covbal2_overall, decimals)),
        adj = 1, line = 3, font = 1)
  
}


fn_vis_stdDiff_wide = function(dat, W=NULL, stra1, # stra2, 
                               ylim=NULL, ylab=NULL, main=NULL, 
                               vertical_ref=FALSE,
                               horizontal_ref=c(-0.1, 0.1),
                               # arrow=TRUE, arrow_length=0.1, 
                               name1=NULL, # name2=NULL,
                               decimals=3 #,  add_legend=FALSE
) {
  # Standardized differences by two subclassifications, no comparison
  # 
  # Parameters
  # ----------
  # W : NULL or character, optional (default=NULL)
  #   The column name of treatment status in dat.
  #   - If W=NULL, it is assumed that dat has a column 'W' for treatment.
  #   - If W is character, it will be used as the name of a column in dat.
  # 
  # Notes
  # -----
  # - This is an alternative version that allow users to use the raw data.
  # [TODO] Allow one stratification proposal
  
  if (is.null(W)) {
    dat_X = dat[, setdiff(colnames(dat), 'W')]
  W = dat[, 'W']
  } else if (is.character(W)) {
    dat_X = dat[, setdiff(colnames(dat), W)]
  W = dat[, W]
  } 
  
  n_covariate = ncol(dat_X)
  
  covbal1 = apply(dat_X, 2, function(one_column){
    return(fn_stdDiff_subAgg(W = W, X_vec = one_column, stra = stra1)$stdDiff_agg)
  })
  # covbal2 = apply(dat_X, 2, function(one_column){
  #   return(fn_stdDiff_subAgg(W = W, X_vec = one_column, stra = stra2)$stdDiff_agg)
  # })
  
  covbal1_overall = fn_stdDiff_overall(covbal1)
  # covbal2_overall = fn_stdDiff_overall(covbal2)
  
  keep_n1 = sum(stra1 != 0)
  # keep_n2 = sum(stra2 != 0)
  discard_n1 = sum(stra1 == 0)
  # discard_n2 = sum(stra2 == 0)
  total_n1 = length(stra1) 
  # total_n2 = length(stra2)
  
  #[TODO] Summary over all covariates
  
  # result_col_list = c('orangered', 'mediumseagreen', 
  #                'lightcoral')
  # result_col = rep('mediumseagreen', length(covbal1))
  # result_col[(abs(covbal2) > abs(covbal1)) & 
  #              (abs(covbal2) <= vertical_ref[2])] = 'lightcoral'
  # result_col[(abs(covbal2) > abs(covbal1)) & 
  #              (abs(covbal2) > vertical_ref[2])] = 'orangered'
  
  if (is.null(ylim)) {
    ylim = max(range(covbal1, na.rm = TRUE)) * c(-1, 1)
    # xlim = max(range(covbal1, covbal2, na.rm = TRUE)) * c(-1, 1)
  }
  
  if (is.null(ylab)) {
    ylab = 'Standardized difference'
  }
  
  if (is.null(main)) {
    main = 'Covariate balance'  
  }
  
  # yticks_at = seq(n_covariate, 1, length.out = 10)
  xticks_at = seq(1, n_covariate, by = 10)
  if (!is.null(names(covbal1))) {
    # yticks = names(covbal1)[seq(1, n_covariate, length.out = 10)]
    xticks = names(covbal1)[seq(1, n_covariate, by = 10)]
  } else {
    # yticks = paste0('X_', seq(1, n_covariate, length.out = 10))
    yticks = paste0('X_', seq(1, n_covariate, by = 10))
  }
  
  plot(y = covbal1, x = 1:n_covariate, pch = 16, # col = 'gray', 
       ylim = ylim, 
       ylab = ylab,
       xaxt = 'n', xlab = '', # 'Covariate', 
       main = main)
  
  axis(1, at = xticks_at, labels = xticks, las = 1)
  abline(h = 0, col = 'gray')
  abline(h = horizontal_ref, col = 'gray', lty = 2)
  
  if (vertical_ref) {
    abline(h = 1:n_covariate, col = 'gray', lty = 2)
  }
  
  # if (add_legend) {
  #   legend('topright', # 'bottomright', 
  #          c('Before', 'Better', 'Acceptable', 'Worse'), 
  #          pch = 16, col = c('gray', 'mediumseagreen', 'lightcoral', 'orangered'))
  # }
  
  title(sub = paste0('Sample sizes: ', keep_n1), adj = 1, line = 2, font = 1)
  title(sub = paste0('Overall balance: ', round(covbal1_overall, decimals)),
        adj = 1, line = 3, font = 1)
  
}