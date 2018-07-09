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

#center variables
ic2surv$NGF_centered <- ic2surv$NGF - mean(ic2surv$NGF)
ic2surv$MAP1B_centered <- ic2surv$MAP1B - mean(ic2surv$MAP1B)

ic2surv <- ic2surv %>%
    dplyr::mutate(mastectomy = breast_surgery == "MASTECTOMY",
                crtherapy = (chemotherapy == "YES" | radio_therapy == "YES"),
                chemotherapy = chemotherapy == "YES",
                radio_therapy = radio_therapy == "YES",
                hormone_therapy = hormone_therapy == "YES",
                kras = ifelse(NGF_centered > 0, T, F),
                lef1 = ifelse(MAP1B_centered > 0, T, F) )
  

#prepare long format dataset
long_ic2dat <- gen_stan_dat(ic2surv)
m1.map2stan <- post_surv(x = ic2surv, surv_form = c("kras", "mastectomy", "kras:mastectomy") )

form <- as.formula( paste0(c( "status~1+offset(log_dtime)+s(time)+", paste(surv_form, collapse = "+")), collapse = ""))

#prepare long format dataset
long_x <- gen_stan_dat(x)

m1.stan_gam <- brms::make_stancode( brms::brm(brms::bf(form),
                 data = long_x, family = poisson() ) )

summary(m1.map2stan)


#Conterfactual plot

timepoints = seq(min(ic2surv$time), max(ic2surv$time), length.out = 100)
dtime = c(0, diff(timepoints))
bcdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  mastectomy = F,
  kras = F
)
bchdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  mastectomy = F,
  kras = T
)
masdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  mastectomy = T,
  kras = F
)
mashdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  mastectomy = T,
  kras = T
)

bc.plot.frame <- get_plot_new_frame(mod = m1.map2stan, newdat = bcdata, treatment = "kras-, underwent breast conservation")
bch.plot.frame <- get_plot_new_frame(mod = m1.map2stan, newdat = bchdata, treatment = "kras+, underwent breast conservation")
mas.plot.frame <- get_plot_new_frame(mod = m1.map2stan, newdat = masdata, treatment = "kras-, underwent mastectomy")
mash.plot.frame <- get_plot_new_frame(mod = m1.map2stan, newdat = mashdata, treatment = "kras+, underwent mastectomy")


data.plot <- rbind(bc.plot.frame, mas.plot.frame, bch.plot.frame, mash.plot.frame)

obs.mortality <- get_km.frame(ic2surv, strata = c("kras", "mastectomy") )
obs.mortality$kras <- as.logical(obs.mortality$kras)
obs.mortality$mastectomy <- as.logical(obs.mortality$mastectomy)

obs.mortality <- obs.mortality %>% 
  dplyr::mutate(Treatment = ifelse(kras & mastectomy, "kras+, underwent mastectomy", ifelse( !(kras) & (mastectomy), "kras-, underwent mastectomy", ifelse( kras &  !mastectomy, "kras+, underwent breast conservation", "kras-, underwent breast conservation"))))

pp <- ggplot2::ggplot(data.plot %>%
                        filter(Treatment == "kras+, underwent mastectomy" | Treatment == "kras+, underwent breast conservation") %>%
                        mutate(Treatment = ifelse(Treatment == "kras+, underwent mastectomy", "mastectomy", "breast conservation") ),
                     aes(time, survmean,
                         group = Treatment, colour = Treatment)) + 
  geom_line(mapping=aes(group = Treatment, colour = Treatment),  alpha = 0.8) +
  geom_ribbon(aes( ymin = survlower,
                   ymax = survupper,
                   fill = Treatment), alpha = 0.5, colour=NA) +
  geom_step(data= obs.mortality %>%
              filter(Treatment == "kras+, underwent mastectomy" | Treatment == "kras+, underwent breast conservation") %>%
              mutate(Treatment = ifelse(Treatment == "kras+, underwent mastectomy", "mastectomy", "breast conservation") ), aes(time,                                                                 surv,                                    colour = Treatment),  alpha = 0.9)+
  labs(title = "iC-2: kras + by type of surgery",
       subtitle = "Survival counterfactual plot",
       caption = "4400 posterior sample size") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Survival probability") +
  xlab("Time (months)") +
  theme_bw() + theme(legend.position=c(0.7, 0.92)) 

pn <- ggplot2::ggplot(data.plot %>%
                        filter(Treatment == "kras-, underwent mastectomy" | Treatment == "kras-, underwent breast conservation") %>%
                        mutate(Treatment = ifelse(Treatment == "kras-, underwent mastectomy", "mastectomy", "breast conservation") ),
                      aes(time, survmean,
                          group = Treatment, colour = Treatment)) + 
  geom_line(mapping=aes(group = Treatment, colour = Treatment),  alpha = 0.8) +
  geom_ribbon(aes( ymin = survlower,
                   ymax = survupper,
                   fill = Treatment), alpha = 0.5, colour=NA) +
  geom_step(data= obs.mortality %>%
              filter(Treatment == "kras-, underwent mastectomy" | Treatment == "kras-, underwent breast conservation") %>%
              mutate(Treatment = ifelse(Treatment == "kras-, underwent mastectomy", "mastectomy", "breast conservation") ), aes(time,                                                                 surv,                                    colour = Treatment),  alpha = 0.9)+
  labs(title = "iC-2: kras - by type of surgery",
       subtitle = "Survival counterfactual plot",
       caption = "4400 posterior sample size") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Survival probability") +
  xlab("Time (months)") +
  theme_bw() + theme(legend.position=c(0.7, 0.92)) 

require(gridExtra)
grid.arrange(pp, pn, ncol = 2)

pdf('counterfactual_kras_surgery.pdf')
grid.arrange(pp, pn, ncol = 2)
dev.off()



#prepare long format dataset
long_ic2dat <- gen_stan_dat(ic2surv)
m2.map2stan <- post_surv(x = ic2surv, surv_form = c("lef1", "mastectomy", "lef1:mastectomy") )

summary(m2.map2stan)


#Conterfactual plot

timepoints = seq(min(ic2surv$time), max(ic2surv$time), length.out = 100)
dtime = c(0, diff(timepoints))
bcdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  mastectomy = F,
  lef1 = F
)
bchdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  mastectomy = F,
  lef1 = T
)
masdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  mastectomy = T,
  lef1 = F
)
mashdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  mastectomy = T,
  lef1 = T
)

bc.plot.frame <- get_plot_new_frame(mod = m2.map2stan, newdat = bcdata, treatment = "lef1-, underwent breast conservation")
bch.plot.frame <- get_plot_new_frame(mod = m2.map2stan, newdat = bchdata, treatment = "lef1+, underwent breast conservation")
mas.plot.frame <- get_plot_new_frame(mod = m2.map2stan, newdat = masdata, treatment = "lef1-, underwent mastectomy")
mash.plot.frame <- get_plot_new_frame(mod = m2.map2stan, newdat = mashdata, treatment = "lef1+, underwent mastectomy")


data.plot <- rbind(bc.plot.frame, mas.plot.frame, bch.plot.frame, mash.plot.frame)

obs.mortality <- get_km.frame(ic2surv, strata = c("lef1", "mastectomy") )
obs.mortality$lef1 <- as.logical(obs.mortality$lef1)
obs.mortality$mastectomy <- as.logical(obs.mortality$mastectomy)

obs.mortality <- obs.mortality %>% 
  dplyr::mutate(Treatment = ifelse(lef1 & mastectomy, "lef1+, underwent mastectomy", ifelse( !(lef1) & (mastectomy), "lef1-, underwent mastectomy", ifelse( lef1 &  !mastectomy, "lef1+, underwent breast conservation", "lef1-, underwent breast conservation"))))

pp <- ggplot2::ggplot(data.plot %>%
                        filter(Treatment == "lef1+, underwent mastectomy" | Treatment == "lef1+, underwent breast conservation") %>%
                        mutate(Treatment = ifelse(Treatment == "lef1+, underwent mastectomy", "mastectomy", "breast conservation") ),
                      aes(time, survmean,
                          group = Treatment, colour = Treatment)) + 
  geom_line(mapping=aes(group = Treatment, colour = Treatment),  alpha = 0.8) +
  geom_ribbon(aes( ymin = survlower,
                   ymax = survupper,
                   fill = Treatment), alpha = 0.5, colour=NA) +
  geom_step(data= obs.mortality %>%
              filter(Treatment == "lef1+, underwent mastectomy" | Treatment == "lef1+, underwent breast conservation") %>%
              mutate(Treatment = ifelse(Treatment == "lef1+, underwent mastectomy", "mastectomy", "breast conservation") ), aes(time,                                                                 surv,                                    colour = Treatment),  alpha = 0.9)+
  labs(title = "iC-2: lef1 + by type of surgery",
       subtitle = "Survival counterfactual plot",
       caption = "4400 posterior sample size") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Survival probability") +
  xlab("Time (months)") +
  theme_bw() + theme(legend.position=c(0.7, 0.92)) 

pn <- ggplot2::ggplot(data.plot %>%
                        filter(Treatment == "lef1-, underwent mastectomy" | Treatment == "lef1-, underwent breast conservation") %>%
                        mutate(Treatment = ifelse(Treatment == "lef1-, underwent mastectomy", "mastectomy", "breast conservation") ),
                      aes(time, survmean,
                          group = Treatment, colour = Treatment)) + 
  geom_line(mapping=aes(group = Treatment, colour = Treatment),  alpha = 0.8) +
  geom_ribbon(aes( ymin = survlower,
                   ymax = survupper,
                   fill = Treatment), alpha = 0.5, colour=NA) +
  geom_step(data= obs.mortality %>%
              filter(Treatment == "lef1-, underwent mastectomy" | Treatment == "lef1-, underwent breast conservation") %>%
              mutate(Treatment = ifelse(Treatment == "lef1-, underwent mastectomy", "mastectomy", "breast conservation") ), aes(time,                                                                 surv,                                    colour = Treatment),  alpha = 0.9)+
  labs(title = "iC-2: lef1 - by type of surgery",
       subtitle = "Survival counterfactual plot",
       caption = "4400 posterior sample size") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Survival probability") +
  xlab("Time (months)") +
  theme_bw() + theme(legend.position=c(0.7, 0.92)) 

require(gridExtra)
grid.arrange(pp, pn, ncol = 2)

pdf('counterfactual_lef1_surgery.pdf')
grid.arrange(pp, pn, ncol = 2)
dev.off()





#prepare long format dataset
long_ic2dat <- gen_stan_dat(ic2surv)
m3.map2stan <- post_surv(x = ic2surv, surv_form = c("kras", "hormone_therapy", "kras:hormone_therapy") )
m4.map2stan <- post_surv(x = ic2surv, surv_form = c("kras", "radio_therapy", "kras:radio_therapy") )
m5.map2stan <- post_surv(x = ic2surv, surv_form = c("kras", "chemotherapy", "kras:chemotherapy") )



summary(m5.map2stan)
summary(m5.map2stan)


#Conterfactual plot

timepoints = seq(min(ic2surv$time), max(ic2surv$time), length.out = 100)
dtime = c(0, diff(timepoints))
bcdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  chemotherapy = F,
  kras = F
)
bchdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  chemotherapy = F,
  kras = T
)
masdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  chemotherapy = T,
  kras = F
)
mashdata <- data.frame(
  time = timepoints,
  log_dtime = log(dtime), 
  chemotherapy = T,
  kras = T
)

bc.plot.frame <- get_plot_new_frame(mod = m5.map2stan, newdat = bcdata, treatment = "kras-, underwent breast conservation")
bch.plot.frame <- get_plot_new_frame(mod = m5.map2stan, newdat = bchdata, treatment = "kras+, underwent breast conservation")
mas.plot.frame <- get_plot_new_frame(mod = m5.map2stan, newdat = masdata, treatment = "kras-, underwent chemotherapy")
mash.plot.frame <- get_plot_new_frame(mod = m5.map2stan, newdat = mashdata, treatment = "kras+, underwent chemotherapy")


data.plot <- rbind(bc.plot.frame, mas.plot.frame, bch.plot.frame, mash.plot.frame)

obs.mortality <- get_km.frame(ic2surv, strata = c("kras", "chemotherapy") )
obs.mortality$kras <- as.logical(obs.mortality$kras)
obs.mortality$chemotherapy <- as.logical(obs.mortality$chemotherapy)

obs.mortality <- obs.mortality %>% 
  dplyr::mutate(Treatment = ifelse(kras & chemotherapy, "kras+, underwent chemotherapy", ifelse( !(kras) & (chemotherapy), "kras-, underwent chemotherapy", ifelse( kras &  !chemotherapy, "kras+, underwent breast conservation", "kras-, underwent breast conservation"))))

pp <- ggplot2::ggplot(data.plot %>%
                        filter(Treatment == "kras+, underwent chemotherapy" | Treatment == "kras+, underwent breast conservation") %>%
                        mutate(Treatment = ifelse(Treatment == "kras+, underwent chemotherapy", "chemotherapy", "no chemotherapy") ),
                      aes(time, survmean,
                          group = Treatment, colour = Treatment)) + 
  geom_line(mapping=aes(group = Treatment, colour = Treatment),  alpha = 0.8) +
  geom_ribbon(aes( ymin = survlower,
                   ymax = survupper,
                   fill = Treatment), alpha = 0.5, colour=NA) +
  geom_step(data= obs.mortality %>%
              filter(Treatment == "kras+, underwent chemotherapy" | Treatment == "kras+, underwent breast conservation") %>%
              mutate(Treatment = ifelse(Treatment == "kras+, underwent chemotherapy", "chemotherapy", "no chemotherapy") ), aes(time,                                                                 surv,                                    colour = Treatment),  alpha = 0.9)+
  labs(title = "iC-2: kras + by type of surgery",
       subtitle = "Survival counterfactual plot",
       caption = "4400 posterior sample size") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Survival probability") +
  xlab("Time (months)") +
  theme_bw() + theme(legend.position=c(0.7, 0.92)) 

pn <- ggplot2::ggplot(data.plot %>%
                        filter(Treatment == "kras-, underwent chemotherapy" | Treatment == "kras-, underwent breast conservation") %>%
                        mutate(Treatment = ifelse(Treatment == "kras-, underwent chemotherapy", "chemotherapy", "no chemotherapy") ),
                      aes(time, survmean,
                          group = Treatment, colour = Treatment)) + 
  geom_line(mapping=aes(group = Treatment, colour = Treatment),  alpha = 0.8) +
  geom_ribbon(aes( ymin = survlower,
                   ymax = survupper,
                   fill = Treatment), alpha = 0.5, colour=NA) +
  geom_step(data= obs.mortality %>%
              filter(Treatment == "kras-, underwent chemotherapy" | Treatment == "kras-, underwent breast conservation") %>%
              mutate(Treatment = ifelse(Treatment == "kras-, underwent chemotherapy", "chemotherapy", "no chemotherapy") ), aes(time,                                                                 surv,                                    colour = Treatment),  alpha = 0.9)+
  labs(title = "iC-2: kras - by type of surgery",
       subtitle = "Survival counterfactual plot",
       caption = "4400 posterior sample size") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Survival probability") +
  xlab("Time (months)") +
  theme_bw() + theme(legend.position=c(0.7, 0.92)) 

require(gridExtra)
grid.arrange(pp, pn, ncol = 2)

pdf('counterfactual_kras_surgery.pdf')
grid.arrange(pp, pn, ncol = 2)
dev.off()
