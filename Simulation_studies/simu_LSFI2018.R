# Lu, Sadiq, Feaster, Ishwaran (2018)
# Estimating individual treatment effect in observational data using random forest methods
# Journal of Computational and Graphical Statistics

# Parameters --------------------------------------------------------------

n_sample = 5000 # 500, 5000
n_simu = 250 # 1000, 250 (corresponding to n_sample)

para_PS_fn = function(X) {
  # input: covariate vector X (length 20)
  plogis( # inverse logit
    -2 + 0.028 * X[1] - 0.374 * X[2] - 0.03 * X[3] + 0.118 * X[4] - 
    0.394 * X[11] + 0.875 * X[12] + 0.9 * X[13]
  ) 
}

para_Y_mean_fn_g = function(X) {
  # input: covariate vector X (length 20)
  0.254 * (X[2]^2) - 0.152 * X[11] - 0.4 * (X[11]^2) - 0.126 * X[12]
}

para_Y_mean_fn_h = function(X) {
  # input: covariate vector X (length 20)
  0.254 * (X[3]^2) - 0.152 * X[4] - 0.126 * X[5] - 0.4 * (X[5]^2)
}

para_Y_mean_fn_1 = function(X, W) {
  # input: covariate vector X (length 20)
  # input: treatment status W (0, 1)
  2.455 + 
    (-1) * (1 - W) * (0.4 * X[1] + 0.154 * X[2] - 0.152 * X[11] - 0.126 * X[12]) + 
    (-1) * W * (para_Y_mean_fn_g(X) > 0) 
}

para_Y_mean_fn_2 = function(X, W) {
  # input: covariate vector X (length 20)
  # input: treatment status W (0, 1)
  2.455 + 
    (-1) * (1 - W) * sin(0.4 * X[1] + 0.154 * X[2] - 0.152 * X[11] - 0.126 * X[12]) + 
    (-1) * W * (para_Y_mean_fn_g(X) > 0)  
}

para_Y_mean_fn_3 = function(X, W) {
  # input: covariate vector X (length 20)
  # input: treatment status W (0, 1)
  2.455 + 
    (-1) * (1 - W) * sin(0.4 * X[1] + 0.154 * X[2] - 0.152 * X[11] - 0.126 * X[12]) + 
    (-1) * W * (para_Y_mean_fn_h(X) > 0)  
}

para_Y_sigma_error = 0.1


# Simulation one dataset --------------------------------------------------

set.seed(20200415 + 2259)

dat_X = cbind(
  matrix(rnorm(n_sample * 11), ncol = 11), 
  matrix(rbinom(n_sample * 9, size = 1, prob = 0.5), ncol = 9)
)
colnames(dat_X) = sprintf("X%02d", 1:20)

dat_PS = apply(dat_X, 1, para_PS_fn)
dat_W = rbinom(n_sample, size = 1, prob = dat_PS)

# Function 1
dat_Y_mean_fn_choice = para_Y_mean_fn_1
dat_Y0 = apply(dat_X, 1, dat_Y_mean_fn_choice, W = 0)
dat_Y1 = apply(dat_X, 1, dat_Y_mean_fn_choice, W = 1)
# dat_Y_obs = mapply(dat_Y_mean_fn_choice, as.data.frame(t(dat_X)), dat_W) + 
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



