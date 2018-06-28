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
m1.map2stan <- post_surv(x = ic2surv, surv_form = c("hormone_therapy", "mastectomy") )

summary(m1.map2stan)


#Conterfactual plot

timepoints = seq(min(ic2surv$time), max(ic2surv$time), length.out = 1000)
dtime = c(0, diff(timepoints))
mastectomy = c( rep(T, 1000), rep(F, 1000))
hormone_therapy = c( rep(F, 1000), rep(T, 1000))
bcdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  mastectomy = F,
  hormone_therapy = T
)
masdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  mastectomy = T,
  hormone_therapy = F
)

bc.plot.frame <- get_plot_new_frame(mod = m1.map2stan, newdat = bcdata, treatment = "Breast conservation with\n hormone therapy")
mas.plot.frame <- get_plot_new_frame(mod = m1.map2stan, newdat = masdata, treatment = "Mastectomy without hormone therapy")


data.plot <- rbind(bc.plot.frame, mas.plot.frame)

obs.mortality <- get_km.frame(ic2surv, strata = c("hormone_therapy", "mastectomy") )

pY <- ggplot2::ggplot(plot.frame ,
                     aes(time, survmean,
                         group = treatment, colour = treatment)) + 
  geom_line(mapping=aes(group = treatment, colour = treatment),  alpha = 0.8) +
  geom_ribbon(aes( ymin = survlower,
                   ymax = survupper,
                   fill = treatment), alpha = 0.5, colour=NA) +
  geom_step(data= obs.mortality %>% 
              dplyr::mutate(treatment = ifelse(as.logical(hormone_therapy) & !as.logical(mastectomy), "Breast conservation with\n hormone therapy", ifelse(!as.logical(hormone_therapy) & as.logical(mastectomy), "Mastectomy without hormone therapy", "something else"))) %>%
              filter(treatment != "something else"), aes(time,                                                                 surv,                                    colour = treatment),  alpha = 0.9)+
  labs(title = "Bayesian Modeling of clinical outcome for iC-2",
       subtitle = "Counterfactual plot by type of surgery and hormone-therapy",
       caption = "4000 posterior samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Survival probability") +
  xlab("Time (months)") +
  theme_bw() + theme(legend.position=c(0.65, 0.85)) 
pdf('counterfactual.pdf')
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

