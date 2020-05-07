# Simulation: AI2016
# Athey, Imbens (2016)
# Recursive partitioning for heterogeneous causal effects
# Proceedings of the National Academy of Sciences of the United States of America


# Overview ----------------------------------------------------------------

# The article proposed the causal tree model. 
# The simulation studies are of small scale, with random treatment 
# assignment (p = 0.5).

set.seed(1)

para_N = 1000


# Designs -----------------------------------------------------------------

#### Design 1 ####

# Mean effect, x = input vector of length K=2
fn_eta1 = function(x) {0.5 * x[1] + x[2]}
# Treatment effect, x = input vector of length K=2
fn_kappa1 = function(x) {0.5 * x[1]}

#### Design 2 ####

# Mean effect, x = input vector of length K=10
fn_eta2 = function(x) {0.5 * sum(x[1:2]) + sum(x[3:6])}
# Treatment effect, x = input vector of length K=10
fn_kappa2 = function(x) {sum((x[1:2] > 0) * x[1:2])}

#### Design 3 ####

# Mean effect, x = input vector of length K=20
fn_eta3 = function(x) {0.5 * sum(x[1:4]) + sum(x[5:8])}
# Treatment effect, x = input vector of length K=20
fn_kappa3 = function(x) {sum((x[1:4] > 0) * x[1:4])}

#### Get X, Y ####

fn_getX = function(size, design) {
  # Get X data 
  if (design == 1) {para_K = 2
  } else if (design == 2) {para_K = 10
  } else if (design == 3) {para_K = 20
  }
  dat_X = matrix(rnorm(size * para_K), ncol = para_K)
  return(dat_X)
}

fn_getY = function(dat_X, dat_W, design) {
  # Get Y data
  if (design == 1) {fn_eta = fn_eta1; fn_kappa = fn_kappa1
  } else if (design == 2) {fn_eta = fn_eta2; fn_kappa = fn_kappa2
  } else if (design == 3) {fn_eta = fn_eta3; fn_kappa = fn_kappa3
  }
  dat_Y = apply(dat_X, 1, fn_eta) + 
          0.5 * (2 * dat_W - 1) * apply(dat_X, 1, fn_kappa) + 
          rnorm(nrow(dat_X), sd = 0.1)
  return(dat_Y)
}

# Get data ----------------------------------------------------------------

prob_W_prob = 0.5
dat_W = rbinom(para_N, size = 1, prob = prob_W_prob)

para_designID = 3
dat_X = fn_getX(size = para_N, design = para_designID)
dat_Y = fn_getY(dat_X = dat_X, dat_W = dat_W, design = para_designID)


# Visualization -----------------------------------------------------------

par(mfrow = c(2, 1))

xLim = range(dat_Y)
hist(dat_Y[dat_W == 0], xlab = 'Y', xlim = xLim, main = 'W == 0')
hist(dat_Y[dat_W == 1], xlab = 'Y', xlim = xLim, main = 'W == 1')

par(mfrow = c(1, 1))
