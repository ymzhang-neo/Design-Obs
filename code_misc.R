# code_misc
# Miscellaneous functions for Design-Obs

# Notes -------------------------------------------------------------------
# 
# Functions in here facilitate the use of other functions in Design-Obs.
# For example, fn_getDummies() will create dummy binary variables for one 
# column of categorical variable.


fn_getDummies = function(x, refLevel = NULL) {
  # Get binary dummy variable for x
  #
  # Parameters
  #   x: a categorical covariate
  #     x could be a factor, a vector of character values, or a vector 
  #     of discrete numeric values.  
  #
  # Returns
  #   result: a matrix of binary dummy variables
  #     The number of columns of the matrix is 1 less than that of the number
  #     of unique values of x.
  #
  # Notes
  # - The function does not check whether x is indeed categorical.
  # - We assume that the original covariate names do not contain string "_dummy_".
  
  if (!is.factor(x)) {
    x = factor(x)
  }
  
  if (!is.null(refLevel)) {
    if (!(refLevel %in% levels(x))) {
      levels(x) = c(refLevel, levels(x))
    } 
    x = relevel(x, refLevel)
  }
  
  result = model.matrix(~ x)[, -1]
  # attr(result, which = 'assign') = NULL
  # attr(result, which = 'contrasts') = NULL
  
  result_colnames = colnames(result)
  colnames(result) = gsub('x', '_dummy_', result_colnames)
  
  return(result)
}