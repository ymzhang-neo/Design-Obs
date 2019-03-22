# learn_MatchIt
# Code for learning to use the package MatchIt (HIKS2011)

library(MatchIt)
data('lalonde')


# Exact matching ----------------------------------------------------------

m_exact = matchit(treat ~ educ + black + hispan, data = lalonde, method = 'exact')
m_exact_data = match.data(m_exact)

# Scratch -----------------------------------------------------------------

table(lalonde$treat)
cbind(sapply(lalonde, class))
demo('exact')
rm(m.out)

str(m_exact_data)
str(lalonde)
table(m_exact_data$subclass)
dim(m_exact_data[m_exact_data$subclass == 1, ]) # 52, 12
table(m_exact_data[m_exact_data$subclass == 1, 'treat']) # N0=13, N1=39

mySubclass = 1
plot(ecdf(m_exact_data$age[(m_exact_data$subclass == mySubclass) & (m_exact_data$treat == 0)]), 
     verticals = TRUE, do.points = FALSE, 
     xlab = 'age', main = paste('Subclass ', mySubclass, 
                                ' (N0=', sum((m_exact_data$subclass == mySubclass) & (m_exact_data$treat == 0)), 
                                ', N1=', sum((m_exact_data$subclass == mySubclass) & (m_exact_data$treat == 1)), 
                                ')', sep = ''))
lines(ecdf(m_exact_data$age[(m_exact_data$subclass == mySubclass) & (m_exact_data$treat == 1)]), 
      verticals = TRUE, do.points = FALSE, col = 'red')
legend('bottomright', c('Control', 'Treated'), col = c('black', 'red'), lwd = 2)
