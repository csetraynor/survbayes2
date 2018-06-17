library(iclust2prog)
library(predsurv)
devtools::document()

data("ic2surv")
?ic2surv

#prepare long format dataset
long_ic2dat <- gen_stan_dat(ic2surv)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())