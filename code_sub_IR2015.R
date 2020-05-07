# code_sub_IR2015
# Functions for subclassification

para_sub_tpval = 0.15
para_sub_size_grp = 3
para_sub_size_str = 3

fn_invLogit = function(x) {
  # Inverse logit transformation
  
  return(1 / (1 + exp(-x)))

}

fn_recTrim = function(dat, colTrt, colV, trim_ctrl_only = FALSE) {
  # Remove the experimental units of one treatment level based on range of the
  # units of the other treatment.
  # 
  # Input:
  #   dat: The data frame under consideration
  #   colTrt: String, column name of the treatment indicator
  #   colV: String, column name of the value for trimming
  #   trim_ctrl_only: logical, 
  #     - trim only control units (and keep all treated units) if TRUE
  #     - trim both control and treated units with extreme colV if FALSE
  # 
  # Output: A list containing two vectors
  #   ind_trim: vector of FALSE, TRUE indicating trimming 
  #             (TRUE for trimming; FALSE for all treated units)
  #   ind_keep: vector of FALSE, TRUE indicating keeping 
  #             (TRUE for keeping; TRUE for all treated units)
  # 
  # Assumptions:
  # - Treatment levels are 0 or 1 (1 for treatment).
  
  ind_trt = dat[, colTrt]
  value_v = dat[, colV]
  
  if (trim_ctrl_only) {
    range_v = range(value_v[ind_trt == 1], na.rm = TRUE)
  } else {
    range_v_c = range(value_v[ind_trt == 0], na.rm = TRUE)
    range_v_t = range(value_v[ind_trt == 1], na.rm = TRUE)
    range_v = c(max(range_v_c[1], range_v_t[1], na.rm = TRUE), min(range_v_c[2], range_v_t[2], na.rm = TRUE))  
  }
  
  ind_keep = ((value_v >= range_v[1]) & (value_v <= range_v[2]))
  ind_trim = !ind_keep
  
  return(list(trim=ind_trim, keep=ind_keep))
  
}

fn_toInt = function(x) {
  # Convert values to integers based on their rank (1 for the smallest).
  # in : a numeric vector
  # out: a numeric vector of integers
  
  x_uniq = unique(x)
  x_uniqrank = rank(x_uniq)
  x_new = x
  for (i in 1:length(x_uniq)) {
    x_new[x == x_uniq[i]] = x_uniqrank[i]
  }
  
  return(x_new)
  
}

fn_check = function(dat, colV, cut_tpval, cut_size_grp, cut_size_str, print_message=TRUE) {
  # For a given stratum, check whether a further split is needed.
  # The function is a modified version of fn_check() in "Code 20181105 stepwise match - simple.R".
  # [TODO] Change dat$W to allow the user-specified treatment covariate.
  
  if ((length(unique(dat[dat$W == 0, colV])) == 1) & 
      (length(unique(dat[dat$W == 1, colV])) == 1)) {
    criterion_t = FALSE
  } else {
    criterion_t = (t.test(dat[, colV] ~ dat$W)$p.value < cut_tpval)
  }
  
  
  Vsplit = median(dat[, colV])
  n_c_l = sum((dat$W == 0) & (dat[, colV] < Vsplit))
  n_t_l = sum((dat$W == 1) & (dat[, colV] < Vsplit))
  n_c_u = sum((dat$W == 0) & (dat[, colV] >= Vsplit))
  n_t_u = sum((dat$W == 1) & (dat[, colV] >= Vsplit))
  
  criterion_size = ((min(n_c_l, n_t_l, n_c_u, n_t_u) > cut_size_grp) & 
                      (min(n_c_l + n_t_l, n_c_u + n_t_u) > cut_size_str))
  
  criterion_split = criterion_t & criterion_size
  
  if (criterion_split) {
    if (print_message) {message("Further split this stratum.")}
    new_Vsplit = Vsplit
    new_stratum = ifelse(dat[, colV] < Vsplit, -0.25, 0.25)
  } else {
    if (print_message) {message("No more split of the current stratum.")}
    new_Vsplit = numeric()
    new_stratum = rep(0, nrow(dat))
  }
  
  return(list(new_Vsplit, new_stratum))
  
}

fn_sub_IR2015 = function(dat, colTrt='W', colPS, 
                         cut_tpval, cut_size_grp, cut_size_str,
                         recTrim=TRUE, trim_ctrl_only = TRUE,
                         print_message=TRUE) {
  # Subclassify data using the procedure in Imben and Rubin (2015)
  
  # Assumptions
  #   - dat contains a column of PS values 
  
  if (recTrim) {
    dat = cbind(dat, 
                stra = as.numeric(
                    fn_recTrim(dat=dat, colTrt=colTrt, colV=colPS, trim_ctrl_only = trim_ctrl_only)$keep
                  )
                )
  } else {
    dat = cbind(dat, stra = 1)
  }
  
  dat_sub = dat[dat$stra == 1, ]
  
  nStra = 1
  PS_cutoff = range(dat_sub[, colPS])
  
  PS_check_split = by(dat_sub, INDICES = dat_sub$stra, 
                      FUN = fn_check, 
                      colV = colPS,
                      cut_tpval = cut_tpval, 
                      cut_size_grp = cut_size_grp,
                      cut_size_str = cut_size_str,
                      print_message = print_message)
                      
  PS_further_split = (sum(sapply(PS_check_split, function(x){length(x[[1]])})) > 0)

  while (PS_further_split) {
  
    # PS_cutoff_new = sapply(PS_check_split, function(x){length(x[[1]])})
  
    for (iStra in 1:nStra) {
      PS_cutoff = c(PS_cutoff, PS_check_split[[iStra]][[1]])
      dat_sub$stra[dat_sub$stra == iStra] = iStra + PS_check_split[[iStra]][[2]]
    }
    dat_sub$stra = fn_toInt(dat_sub$stra)
    
    nStra = length(PS_cutoff) - 1
    
    PS_check_split = by(dat_sub, INDICES = dat_sub$stra, 
                        FUN = fn_check, 
                        colV = colPS,
                        cut_tpval = cut_tpval, 
                        cut_size_grp = cut_size_grp,
                        cut_size_str = cut_size_str,
                        print_message = print_message)
    
    PS_further_split = (sum(sapply(PS_check_split, function(x){length(x[[1]])})) > 0)
  
  }
  
  dat$stra[dat$stra == 1] = dat_sub$stra
  
  return(list(stra=dat$stra, cutoff=PS_cutoff))
  
}



