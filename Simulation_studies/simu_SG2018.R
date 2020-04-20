# Samuels, Greevy (2018)
# Bagged one-to-one matching for efficient and robust treatment effect estimation
# Statistics in Medicine

# Parameters --------------------------------------------------------------

n_sample = 1000

para_X_p = 10

para_PS_W1_target = 0.1
para_beta0 = log(para_PS_W1_target / (1 - para_PS_W1_target))
# para_beta0 = -3.05
para_beta_X = c(0.25, 0.50, 0.75, 0.90)
names(para_beta_X) = c('L', 'M', 'H', 'VH')

para_PS_fn = function(X) {
  # input: covariate vector X 
  plogis(
    para_beta0 + 
      para_beta_X['L'] * X[1] + 
      para_beta_X['M'] * X[2] + 
      para_beta_X['H'] * X[3] + 
      0.5 * para_beta_X['L'] * X[4] + 
      0.3 * para_beta_X['L'] * (X[4]^2) + 
      0.5 * para_beta_X['M'] * X[5] + 
      0.5 * para_beta_X['H'] * X[6] + 
      0.3 * para_beta_X['H'] * X[5] * X[6] + 
      0.5 * para_beta_X['VH'] * X[7] + 
      0.3 * para_beta_X['VH'] * (X[7]^2)
  )
}

para_Y_TE = 1
para_Y_mean_fn = function(X, W) {
  # input: covariate vector X (length 20)
  # input: treatment status W (0, 1)
  para_Y_TE * W + 
    0.5 * para_beta_X['L'] * X[4] + 
    0.3 * para_beta_X['L'] * (X[4]^2) + 
    0.5 * para_beta_X['M'] * X[5] + 
    0.5 * para_beta_X['H'] * X[6] + 
    0.3 * para_beta_X['H'] * X[5] * X[6] + 
    0.5 * para_beta_X['VH'] * X[7] + 
    0.3 * para_beta_X['VH'] * (X[7]^2) + 
    para_beta_X['L'] * X[8] + 
    para_beta_X['M'] * X[9] + 
    para_beta_X['H'] * X[10]
}
para_Y_sigma_error = 3


# Simulate one dataset ----------------------------------------------------

set.seed(20200416 + 0111)

dat_X = matrix(rnorm(n_sample * para_X_p), ncol = para_X_p)
colnames(dat_X) = sprintf("X%02d", 1:para_X_p)

dat_PS = apply(dat_X, 1, para_PS_fn)
dat_W = rbinom(n_sample, size = 1, prob = dat_PS)

# Function 1
dat_Y0 = apply(dat_X, 1, para_Y_mean_fn, W = 0)
dat_Y1 = apply(dat_X, 1, para_Y_mean_fn, W = 1)
# dat_Y_obs = mapply(para_Y_mean_fn, as.data.frame(t(dat_X)), dat_W) + 
#   rnorm(n_sample, sd = para_Y_sigma_error)

dat_Y_obs = (1 - dat_W) * dat_Y0 + dat_W * dat_Y1 + 
  rnorm(n_sample, sd = para_Y_sigma_error)
dat_Y = cbind(dat_Y0, dat_Y1, dat_Y_obs)
colnames(dat_Y) = c('Y0', 'Y1', 'Y')

dat = data.frame(dat_X, PS = dat_PS, W = dat_W, dat_Y)


# Estimate without design (D0) --------------------------------------------

D00_est_ATE = mean(dat[, 'Y1']) - mean(dat[, 'Y0'])

#### D0 Neyman ####

D0_estN_est = mean(dat[dat$W == 1, 'Y']) - mean(dat[dat$W == 0, 'Y'])

D0_estN_var_W1 = var(dat[dat$W == 1, 'Y'])
D0_estN_var_W0 = var(dat[dat$W == 0, 'Y'])
D0_estN_var = D0_estN_var_W1 / sum(dat$W == 1) + 
  D0_estN_var_W0 / sum(dat$W == 0)

D0_estN_CI = D0_estN_est + c(-1, 1) * 1.96 * sqrt(D0_estN_var)

#### D0 Regression (linear, main) ####

D0_estR1_fit = lm(Y ~ . - PS - Y0 - Y1, data = dat)
D0_estR1_est = coef(summary(D0_estR1_fit))['W', 'Estimate']
D0_estR1_CI = D0_estR1_est + c(-1, 1) * 
  1.96 * coef(summary(D0_estR1_fit))['W', 'Std. Error']


# Notes -------------------------------------------------------------------

# References:
# - AS2014
# - Supplement 2 R code: https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fsim.7926&file=SIM_7926-SUP-0002.pdf
# - GitHub code: https://github.com/LaurenSamuels/BOOM