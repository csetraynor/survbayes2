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


m1.gam <- mgcv::gam(status ~ 1 + offset(log_dtime) + s(time) + factor(chemotherapy)*factor(radio_therapy), long_ic2dat, family='poisson')

m1.stan_gam <- rstanarm::stan_gamm4(status~1+offset(log_dtime)+s(time) + factor(chemotherapy), data = long_ic2dat , family='poisson')
summary(m1.stan_gam)

#plot predictions
post <- posterior_linpred(m1.stan_gam)

plot.matrix <- link.surv(post = post, longdata = long_ic2dat)
surv.mean <- as.data.frame( t( apply(plot.matrix, 1, mean) ), colnames = surv.mean)
surv.pi <- as.data.frame( t( apply(plot.matrix, 1, quantile, probs = c(0.055, 0.945)) ), colnames = c( "lower", "upper") )

plot.frame <- get_plot.frame(post = plot.matrix, strata = c("chemotherapy", "radio_therapy") )
obs.mortality <- get_km.frame(ic2surv, strata = c("chemotherapy", "radio_therapy") )

p <- ggplot2::ggplot(plot.frame ,
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
m2.stan_gam <- rstanarm::stan_gamm4(status ~ 1 + offset(log_dtime) + s(time) + factor(chemotherapy)*factor(radio_therapy) , data = long_ic2dat , family='poisson')
summary(m2.stan_gam)

post <- posterior_linpred(m2.stan_gam)

plot.matrix <- link.surv(post = post, longdata = long_ic2dat)
surv.mean <- as.data.frame( t( apply(plot.matrix, 1, mean) ), colnames = surv.mean)
surv.pi <- as.data.frame( t( apply(plot.matrix, 1, quantile, probs = c(0.055, 0.945)) ), colnames = c( "lower", "upper") )

plot.frame <- get_plot.frame(post = plot.matrix, strata = c("chemotherapy", "radio_therapy") )
obs.mortality <- get_km.frame(ic2surv, strata = c("chemotherapy", "radio_therapy") )


pNO <- ggplot2::ggplot(plot.frame %>% 
                       filter(radio_therapy == "NO"),
                     aes(time, postmean, group = chemotherapy, colour = chemotherapy)) + 
  geom_line(mapping=aes(group = chemotherapy, colour = chemotherapy),  alpha = 0.8) +
  geom_ribbon(aes( ymin = lower,
                   ymax = upper,
                   fill = chemotherapy), alpha = 0.5, colour=NA) +
  geom_step(data=obs.mortality %>%
              filter(radio_therapy == "NO"), aes(time,
                                                                                                    surv,
                                                                                                    colour = chemotherapy),  alpha = 0.9)+
  ggtitle("No radiotherapy") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Survival probability") +
  xlab("Time (months)") +
  theme_bw() + theme(legend.position=c(0.75, 0.75))

pYES <- ggplot2::ggplot(plot.frame %>% 
                         filter(radio_therapy == "YES"),
                       aes(time, postmean, group = chemotherapy, colour = chemotherapy)) + 
  geom_line(mapping=aes(group = chemotherapy, colour = chemotherapy),  alpha = 0.8) +
  geom_ribbon(aes( ymin = lower,
                   ymax = upper,
                   fill = chemotherapy), alpha = 0.5, colour=NA) +
  geom_step(data=obs.mortality %>%
              filter(radio_therapy == "YES"), aes(time,
                                                 surv,
                                                 colour = chemotherapy),  alpha = 0.9)+
  ggtitle("Yes radiotherapy") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Survival probability") +
  xlab("Time (months)") +
  theme_bw() + theme(legend.position=c(0.75, 0.75))

require(gridExtra)

pdf('TreatEffectInteraction.pdf')
grid.arrange(pNO, pYES, ncol=2)
dev.off()

