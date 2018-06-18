library(iclust2prog)
library(predsurv)
library(magrittr)
library(dplyr)
library(forcats)
library(tidyr)
library(purrr)
library(modelr)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(ggridges)
library(rstan)
library(rstanarm)
import::from(LaplacesDemon, invlogit)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
devtools::document()

data("ic2surv")
?ic2surv

#prepare long format dataset
long_ic2dat <- gen_stan_dat(ic2surv)

m1.stan_gam <- rstanarm::stan_gamm4(status~1+offset(log_dtime)+s(time) + factor(chemotherapy) , data = long_ic2dat , family='poisson')
summary(m1.stan_gam)

#plot predictions
post <- posterior_predict(m1.stan_gam)

plot.list <- lapply(1:nrow(post), function(i){
  long_ic2dat$test <- post[i, ]
  p1 <- predict(mgcv::gam(test ~ 1 + offset(log_dtime) + s(time) + factor(chemotherapy), long_ic2dat, family='poisson'), long_ic2dat)
  long_ic2dat$loghaz <- p1
  dat <- long_ic2dat %>%
    group_by(sample_id) %>%
    mutate(surv = exp( - cumsum(exp(loghaz))) ) %>%
    filter(row_number()==n()) %>%
    ungroup() %>%
    arrange(time) %>%
    select(surv)
  return(dat)
})
plot.matrix <- do.call(cbind, plot.list)
plot.matrix <- as.matrix(plot.matrix)

surv.mean <- as.data.frame( t( apply(plot.matrix, 1, mean) ), colnames = surv.mean)
surv.pi <- as.data.frame( t( apply(plot.matrix, 1, quantile, probs = c(0.055, 0.945)) ), colnames = c( "lower", "upper") )

plot.frame <- data.frame(
  postmean = apply(plot.matrix, 1, mean),
  lower = apply(plot.matrix, 1, quantile, probs = 0.055),
  upper = apply(plot.matrix, 1, quantile, probs = 0.945),
  time = ic2surv %>% arrange(time) %>% select(time) %>% 
                           unlist %>% as.double,
  strata = ic2surv %>% arrange(time) %>% select(chemotherapy) %>% 
    unlist %>% as.character()
  )

mle.surv <- survival::survfit(Surv( time , status) ~ chemotherapy,
                              data = ic2surv  )
obs.mortality <- data.frame(time = mle.surv$time,
                            surv = mle.surv$surv,
                            strata = summary(mle.surv, censored=T)$strata)
ggplot2::ggplot(plot.frame, aes(time, postmean, col = "strata" )) + 
  geom_line() +
  geom_ribbon(aes( ymin = lower,
                   ymax = upper,
                   group = strata), alpha = 0.2) +
  geom_step(data = obs.mortality, aes(time,
                                      surv,
                                      group = strata) )

#Conterfactual plot
m2.stan_gam <- rstanarm::stan_gamm4(status~1+offset(log_dtime)+s(time)+factor(chemotherapy)+factor(radio_therapy)+age_std+npi+NGF+MAP1B, data = long_ic2dat , family='poisson')
summary(m1.stan_gam)
