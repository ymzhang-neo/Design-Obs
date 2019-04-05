# learn_MatchIt
# Code for learning to use the package MatchIt (HIKS2011)

library(MatchIt)
data('lalonde')


# Exact matching ----------------------------------------------------------

m_exact = matchit(treat ~ educ + black + hispan, 
                  data = lalonde, method = 'exact')
m_exact_data = match.data(m_exact)

str(m_exact_data) # 560, 12
str(lalonde) # 614, 10
table(m_exact_data$subclass)

mySubclass = 1

dim(m_exact_data[m_exact_data$subclass == mySubclass, ]) # 52, 12
table(m_exact_data[m_exact_data$subclass == mySubclass, 'treat']) # N0=13, N1=39

plot(ecdf(m_exact_data$age[(m_exact_data$subclass == mySubclass) & (m_exact_data$treat == 0)]), 
     verticals = TRUE, do.points = FALSE, 
     xlab = 'age', main = paste('Subclass ', mySubclass, 
                                ' (N0=', sum((m_exact_data$subclass == mySubclass) & (m_exact_data$treat == 0)), 
                                ', N1=', sum((m_exact_data$subclass == mySubclass) & (m_exact_data$treat == 1)), 
                                ')', sep = ''))
lines(ecdf(m_exact_data$age[(m_exact_data$subclass == mySubclass) & (m_exact_data$treat == 1)]), 
      verticals = TRUE, do.points = FALSE, col = 'red')
legend('bottomright', c('Control', 'Treated'), col = c('black', 'red'), lwd = 2)


# Nearest neighbor matching -----------------------------------------------

m_nn = matchit(treat ~ re74 + re75 + educ + black + hispan + age, 
               data = lalonde, method = 'nearest')
m_nn_data = match.data(m_nn)

str(m_nn_data) # 370, 12


# Optimal matching --------------------------------------------------------

m_opt = matchit(treat ~ re74 + re75 + age + educ,
                data = lalonde, method = 'optimal', ratio = 2)
# Warning message:
# In optmatch::fullmatch(d, min.controls = ratio, max.controls = ratio,  :
#   Without 'data' argument the order of the match is not guaranteed
#     to be the same as your original data.


# Scratch -----------------------------------------------------------------

set.seed(20190322.0037)
tempVec = rbinom(n = 3, size = 1, prob = 0.5)
tempMat = matrix(rbinom(n = 3*40, size = 1, prob = 0.5), ncol = 3)

all.equal(tempMat, tempVec)
apply(tempMat, 1, all.equal, current = tempVec)
# Try this with the actual data