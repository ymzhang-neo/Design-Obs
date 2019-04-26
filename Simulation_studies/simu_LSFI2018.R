# Simulation: LSFI2018
# Lu, Sadiq, Feaster, Ishwaran (2018)
# Estimating individual treatment effect in observational data using random forest methods
# Journal of Computational and Graphical Statistics

set.seed(1)

para_N = 500 # 5000

para_beta_W = c(-2, 0.028, -0.374, -0.03, 0.118, -0.394, 0.875, 0.9)

dat_X = cbind(
  matrix(rnorm(para_N * 11), ncol = 11), 
  matrix(rbinom(para_N * 9, size = 1, prob = 0.5), ncol = 9)
)
colnames(dat_X) = paste0('X', 1:20)

dat_W = rbinom(para_N, size = 1, 
               prob = 1/(1 + exp(- (dat_X[, paste0('X', c(1, 2, 3, 4, 11, 12, 13))] %*% para_beta_W[-1] + para_beta_W[1]))))


# Simulate outcomes -------------------------------------------------------

para_epsilon_sd = 0.1

fn_Y1 = function(X = dat_X, W = dat_W) {
  myY = rep(NA, para_N)
  myY[dat_W == 0] = rnorm(para_N, sd = para_epsilon_sd) + 
    2.455 - dat_X[, paste0('X', c(1, 2, 11, 12))] %*% c(0.4, 0.154, -0.152, -0.126)
  myY[dat_W == 1] = rnorm(para_N, sd = para_epsilon_sd) + 
    2.455 - ((0.254*(dat_X[, 'X2'])^2 - 0.152*dat_X[, 'X11'] - 0.4*dat_X[, 'X11']^2 - 0.126*dat_X[, 'X12']) > 0)
  return(myY)
}

fn_Y2 = function(X = dat_X, W = dat_W) {
  myY = rep(NA, para_N)
  myY[dat_W == 0] = rnorm(para_N, sd = para_epsilon_sd) + 
    2.455 - sin(dat_X[, paste0('X', c(1, 2, 11, 12))] %*% c(0.4, 0.154, -0.152, -0.126))
  myY[dat_W == 1] = rnorm(para_N, sd = para_epsilon_sd) + 
    2.455 - ((0.254*(dat_X[, 'X2'])^2 - 0.152*dat_X[, 'X11'] - 0.4*dat_X[, 'X11']^2 - 0.126*dat_X[, 'X12']) > 0)
  return(myY)
}

fn_Y3 = function(X = dat_X, W = dat_W) {
  myY = rep(NA, para_N)
  myY[dat_W == 0] = rnorm(para_N, sd = para_epsilon_sd) + 
    2.455 - sin(dat_X[, paste0('X', c(1, 2, 11, 12))] %*% c(0.4, 0.154, -0.152, -0.126))
  myY[dat_W == 1] = rnorm(para_N, sd = para_epsilon_sd) + 
    2.455 - ((0.254*(dat_X[, 'X3'])^2 - 0.152*dat_X[, 'X4'] - 0.4*dat_X[, 'X5']^2 - 0.126*dat_X[, 'X5']) > 0)
  return(myY)
}

dat_Y = fn_Y1(X = dat_X, W = dat_W)
# dat_Y = fn_Y2(X = dat_X, W = dat_W)
# dat_Y = fn_Y3(X = dat_X, W = dat_W)