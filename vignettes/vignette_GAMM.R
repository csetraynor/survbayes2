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
library(brms)
library(brmstools)
import::from(LaplacesDemon, invlogit)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
devtools::document()

data("ic2surv")
?ic2surv

ic2surv <- ic2surv %>%
    dplyr::mutate(mastectomy = breast_surgery == "MASTECTOMY",
                crtherapy = (chemotherapy == "YES" | radio_therapy == "YES"),
                chemotherapy = chemotherapy == "YES",
                radio_therapy = radio_therapy == "YES",
                hormone_therapy = hormone_therapy == "YES")
  

#prepare long format dataset
long_ic2dat <- gen_stan_dat(ic2surv)


m1.gam <- mgcv::gam(status ~ 1 + offset(log_dtime) + s(time), long_ic2dat, family='poisson')

m1.stan_gam <- rstanarm::stan_gamm4(status~1+offset(log_dtime)+s(time) + chemotherapy, data = long_ic2dat, family='poisson')
prior <- get_prior(m1.stan_gam)

m1.brms_gam <- brms::brm(bf(status ~ 1 + offset(log_dtime) + s(time) + crtherapy),
                   data = long_ic2dat, family = poisson(), cores = 4, seed = 17,
                   iter = 4000, warmup = 1000, thin = 10, refresh = 0,
                   control = list(adapt_delta = 0.99))


the_stan_model <- make_stancode(bf(status ~ 1 + offset(log_dtime) + s(time) + crtherapy),
                         data = long_ic2dat, family = poisson())



(prior <- get_prior(brm(bf(status ~ 1 + offset(log_dtime) + s(time) + factor(chemotherapy)),
                        data = long_ic2dat, family = poisson())))

summary(m1.brms_gam)

#plot predictions
timepoints <-  seq(0, max(ic2surv$time),
                   length.out = 100L)

newdata <- gen_new.frame(ic2surv, timepoints = timepoints)

chemdata <- rbind(newdata %>% mutate(chemotherapy = T) ,                     newdata %>% mutate(chemotherapy = F,
                  sample_id = sample_id+nrow(ic2surv))                                                             ) %>%
  select(-status)
post <- posterior_linpred(m1.stan_gam, newdata = chemdata)

plot.matrix <- link.surv(post = post, longdata = chemdata)
surv.mean <- as.data.frame( t( apply(plot.matrix, 1, mean) ), colnames = surv.mean)
surv.pi <- as.data.frame( t( apply(plot.matrix, 1, quantile, probs = c(0.055, 0.945)) ), colnames = c( "lower", "upper") )

plot.frame <- get_plot.frame(post = plot.matrix, strata = c("chemotherapy"), obs = rbind(ic2surv, ic2surv))
obs.mortality <- get_km.frame(ic2surv, strata = c("chemotherapy") )

p <- ggplot2::ggplot(plot.frame ,
                     aes(time, postmean, group = chemotherapy, colour = chemotherapy)) + 
  geom_line(mapping=aes(group = chemotherapy, colour = chemotherapy),  alpha = 0.8) +
  geom_ribbon(aes( ymin = lower,
                   ymax = upper,
                   fill = chemotherapy), alpha = 0.5, colour=NA) +
  geom_step(data=obs.mortality , aes(time,
                                    surv,
                                    colour = chemotherapy),  alpha = 0.9)+
  labs(title="Bayesian Modeling of survival outcome for iC-2",
       subtitle="Hormone therapy observed survival and predictive posterior distribution with 89% credible interval",
       caption="4000 posterior sample smooths from GAMM model shown") +
  ylab("Survival probability") +
  xlab("Time (months)") +
  # guides() +
  theme_bw() + theme(legend.position=c(0.75, 0.75))
p
pdf('Hormone.pdf')
p
dev.off()





#Conterfactual plot
ic2surv <- ic2surv %>%
  mutate(Treatment = ifelse(chemotherapy & !hormone_therapy,
                "only_chemo", ifelse(hormone_therapy & !chemotherapy, "only_hormone", ifelse(hormone_therapy & chemotherapy, "both", "none"))))


m2.stan_gam <- rstanarm::stan_gamm4(status~1+offset(log_dtime)+s(time) + mastectomy + hormone_therapy, data = long_ic2dat, family='poisson', seed = 17,iter = 4000, warmup = 1000, thin = 10, control = list(adapt_delta = 0.99))


summary(m2.stan_gam)




post <- posterior_linpred(m2.stan_gam)

plot.matrix <- link.surv(post = post, longdata = long_ic2dat)
surv.mean <- as.data.frame( t( apply(plot.matrix, 1, mean) ), colnames = surv.mean)
surv.pi <- as.data.frame( t( apply(plot.matrix, 1, quantile, probs = c(0.055, 0.945)) ), colnames = c( "lower", "upper") )

plot.frame <- get_plot.frame(post = plot.matrix, strata = c("hormone_therapy", "mastectomy") ) 
obs.mortality <- get_km.frame(ic2surv, strata = c("hormone_therapy", "mastectomy") )




pY <- ggplot2::ggplot(plot.frame %>% 
                       dplyr::mutate(Treatment = ifelse(hormone_therapy & !mastectomy, "Breast conservation with\n hormone therapy", ifelse(!hormone_therapy & mastectomy, "Mastectomy without hormone therapy", "something else"))) %>%
                        filter(Treatment != "something else"),
                     aes(time, postmean, group = Treatment, colour = Treatment)) + 
  geom_line(mapping=aes(group = Treatment, colour = Treatment),  alpha = 0.8) +
  geom_ribbon(aes( ymin = lower,
                   ymax = upper,
                   fill = Treatment), alpha = 0.5, colour=NA) +
  geom_step(data= obs.mortality %>% 
              dplyr::mutate(Treatment = ifelse(as.logical(hormone_therapy) & !as.logical(mastectomy), "Breast conservation with\n hormone therapy", ifelse(!as.logical(hormone_therapy) & as.logical(mastectomy), "Mastectomy without hormone therapy", "something else"))) %>%
              filter(Treatment != "something else"), aes(time,
                                                                                                    surv,
                                                                                                    colour = Treatment),  alpha = 0.9)+
  labs(title="Bayesian Modeling of clinical outcome for iC-2",
       subtitle="Observed vs predicted posterior survival distribution with 89% credible interval\n by type of surgery and hormone therapy") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Survival probability") +
  xlab("Time (months)") +
  theme_bw() + theme(legend.position=c(0.65, 0.85))
pdf('HormoneSurgeryInteraction.pdf')
pY
dev.off()


pN <- ggplot2::ggplot(plot.frame %>% 
                        dplyr::filter(!mastectomy),
                      aes(time, postmean, group = hormone_therapy, colour = hormone_therapy)) + 
  geom_line(mapping=aes(group = hormone_therapy, colour = hormone_therapy),  alpha = 0.8) +
  geom_ribbon(aes( ymin = lower,
                   ymax = upper,
                   fill = hormone_therapy), alpha = 0.5, colour=NA) +
  geom_step(data=obs.mortality %>%
              
              dplyr::filter(!as.logical(mastectomy)), aes(time,
                                                         surv,
                                                         colour = hormone_therapy),  alpha = 0.9)+
  ggtitle("Practised breast conservation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Survival probability") +
  xlab("Time (months)") +
  theme_bw() + theme(legend.position=c(0.75, 0.75))
pN
require(gridExtra)
pdf('HormoneRadioInteraction.pdf')
grid.arrange(pN, pY, ncol=2)
dev.off()

