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


m1.gam <- mgcv::gam(status ~ 1 + offset(log_dtime) + s(time) + factor(chemotherapy), long_ic2dat, family='poisson')

m1.stan_gam <- rstanarm::stan_gamm4(status~1+offset(log_dtime)+s(time) + factor(chemotherapy) , data = long_ic2dat , family='poisson')
summary(m1.stan_gam)

#plot predictions
post <- posterior_linpred(m1.stan_gam)

plot.list <- lapply(1:nrow(post), function(i){
  long_ic2dat$loghaz <- post[i, ]
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


# surv.mean <- as.data.frame( t( apply(plot.matrix, 1, mean) ), colnames = surv.mean)
# surv.pi <- as.data.frame( t( apply(plot.matrix, 1, quantile, probs = c(0.055, 0.945)) ), colnames = c( "lower", "upper") )

plot.frame <- data.frame(
  postmean = apply(plot.matrix, 1, mean),
  lower = apply(plot.matrix, 1, quantile, probs = 0.055),
  upper = apply(plot.matrix, 1, quantile, probs = 0.945),
  time = ic2surv %>% arrange(time) %>% select(time) %>% 
                           unlist %>% as.double,
  strata = ic2surv %>% arrange(time) %>% select(chemotherapy) %>% 
    unlist %>% as.character()
  )

zeros.frame <- data.frame(time=0, postmean=1, strata=unique(plot.frame$strata),
                          lower = 1, upper = 1)
plot.frame <- rbind(plot.frame, zeros.frame)
mle.surv <- survival::survfit(Surv( time , status) ~ chemotherapy,
                              data = ic2surv  )
obs.mortality <- data.frame(time = mle.surv$time,
                            surv = mle.surv$surv,
                            strata = summary(mle.surv, censored=T)$strata)
zeros <- data.frame(time=0, surv=1, strata=unique(obs.mortality$strata))
obs.mortality <- rbind(obs.mortality, zeros)




p <- ggplot2::ggplot(plot.frame %>% 
                       dplyr::mutate(chemotherapy = strata),
                     aes(time, postmean, group = chemotherapy, colour = chemotherapy)) + 
  geom_line(mapping=aes(group = chemotherapy, colour = chemotherapy),  alpha = 0.8) +
  geom_ribbon(aes( ymin = lower,
                   ymax = upper,
                   fill = chemotherapy), alpha = 0.5, colour=NA) +
  geom_step(data=obs.mortality %>% 
              dplyr::mutate(chemotherapy = ifelse(strata == "chemotherapy=YES", "YES", "NO") ), aes(time,
                                    surv,
                                    colour = chemotherapy),  alpha = 0.9)+
  labs(title="Bayesian Modeling of survival outcome for iC-2",
       subtitle="Comparison of observed survival and predictive posterior distribution with 89% credible interval ",
       caption="4000 posterior sample smooths from GAMM model shown") +
  ylab("Survival probability") +
  xlab("Time (months)") +
  # guides() +
  theme_bw() + theme(legend.position=c(0.75, 0.75))

pdf('survbayesGAMM.pdf')
p
dev.off()





#Conterfactual plot
m2.stan_gam <- rstanarm::stan_gamm4(status~1+offset(log_dtime)+s(time)+factor(chemotherapy)+factor(radio_therapy)+age_std+npi+NGF+MAP1B, data = long_ic2dat , family='poisson')
summary(m1.stan_gam)

post <- as.matrix(m1.stan_gam)
pairs(m1.stan_gam, regex_pars = "^s")
plot(m1.stan_gam, regex_pars = "^s")



