pacman::p_load(dplyr, readr)
library(MiaoCom)

d = read_csv('misc/dataTest.csv', locale=locale(encoding = 'UTF8')) %>%
  select(12:21, age) %>%
  `colnames<-`(c(paste0('diag', 1:10), 'age'))
d1 = cci(d, paste0('diag', 1:10), 'age')
d2 = eci(d, paste0('diag', 1:10))
