# Simulation: WA2018
# Wager, Athey (2018)
# Estimation and inference of heterogeneous treatment effects using random forests
# Journal of the American Statistical Association

# Overview ----------------------------------------------------------------



#### Design 1 ####

para_d1 = 2 # 5, 10, 15, 20, 30

para_N1 = 500 

# PS model, x = input variable (X1)
fn_e1 = function(x) {(1 + dbeta(x, shape1 = 2, shape2 = 4)) / 4}
# Main effect, x = input variable (X1)
fn_m1 = function(x) {2 * x - 1}

fn_tau1 = function(x) {0}

#### Design 2 ####

para_d2 = 2 # 3, 4, 5, 6, 8

para_N2 = 5000

fn_e2 = function(x) {0.5}
fn_m2 = function(x) {0}

# Treatment effect, x = input vector of length d
fn_tau2 = function(x) {
  (1 + 1 / (1 + exp(-20 * (x[1] - 1 / 3)))) * 
  (1 + 1 / (1 + exp(-20 * (x[2] - 1 / 3))))
} 

#### Design 3 ####

para_d3 = 2 # 3, 4, 5, 6, 8

para_N3 = 10000

fn_e3 = function(x) {0.5}
fn_m3 = function(x) {0}

# Treatment effect, x = input vector of length d
fn_tau3 = function(x) {
  (1 + 2 / (1 + exp(-12 * (x[1] - 1 / 2)))) * 
  (1 + 2 / (1 + exp(-12 * (x[2] - 1 / 2))))
} 

#### Design 4 #### (in appendix)

para_d4 = 6 # 12

para_N4 = 5000

fn_e4 = function(x) {0.5}
fn_m4 = function(x) {0}

# Treatment effect, x = input vector of length d
# q = number of covariates with signal, 2, 4, 6
fn_varsigma4 = function(x) {(1 / (1 + exp(-12 * (x - 0.5))) - 0.5)}
fn_tau4 = function(x, q) {sum(sapply(x[1:q], fn_varsigma4)) * 4 / q}

#### Get X, W, Y ####

fn_getX = function(design) {
  # Get X data (design = 1, 2, 3, 4)
  if (design == 1) {size = para_N1; para_d = para_d1
  } else if (design == 2) {size = para_N2; para_d = para_d2
  } else if (design == 3) {size = para_N3; para_d = para_d3
  } else if (design == 4) {size = para_N4; para_d = para_d4
  }
  dat_X = matrix(runif(size * para_d), ncol = para_d)
  return(dat_X)
}

fn_getW = function(dat_X, design) {
  if (design == 1) {fn_e = fn_e1
  } else if (design == 2) {fn_e = fn_e2
  } else if (design == 3) {fn_e = fn_e3
  } else if (design == 4) {fn_e = fn_e4
  }
  dat_W = rbinom(n = nrow(dat_X), size = 1, prob = apply(dat_X, 1, fn_e))
  return(dat_W)
}

fn_getY = function(dat_X, dat_W, design) {
  if (design == 1) {fn_e = fn_e1; fn_m = fn_m1; fn_tau = fn_tau1
  } else if (design == 2) {fn_e = fn_e2; fn_m = fn_m2; fn_tau = fn_tau2
  } else if (design == 3) {fn_e = fn_e3; fn_m = fn_m3; fn_tau = fn_tau3
  } else if (design == 4) {fn_e = fn_e4; fn_m = fn_m4; fn_tau = fn_tau4
  }
  dat_Y = apply(dat_X, 1, fn_m) + 
          0.5 * (2 * dat_W - 1) * apply(dat_X, 1, fn_tau) + 
          rnorm(nrow(dat_X), sd = 1)
}


# Get data ----------------------------------------------------------------

para_designID = 1
dat_X = fn_getX(design = para_designID)
dat_W = fn_getW(dat_X = dat_X, design = para_designID)
dat_Y = fn_getY(dat_X = dat_X, dat_W = dat_W, design = para_designID)
# Notes
# - The current code still requires the users to manually specify
#   para_d1, para_d2, para_d3, para_d4, para_N2.
# - The current code does not allow easy specification of q in Design 4.


# Visualization -----------------------------------------------------------

par(mfrow = c(2, 1))

xLim = range(dat_Y)
hist(dat_Y[dat_W == 0], xlab = 'Y', xlim = xLim, main = 'W == 0')
hist(dat_Y[dat_W == 1], xlab = 'Y', xlim = xLim, main = 'W == 1')

par(mfrow = c(1, 1))