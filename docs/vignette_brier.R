library(iclust2prog)
library(predsurv)
devtools::document()

# Calculta Brier Score ----------------
long_ic2dat <- gen_stan_dat(ic2surv)
set.seed(10)
m2.stan_gam <- rstanarm::stan_gamm4(status ~ 1 + offset(log_dtime) + s(time) + factor(chemotherapy)*factor(radio_therapy) , data = long_ic2dat  , family='poisson')


timepoints <-  seq(0, max(ic2surv$time),
                   length.out = 100L)

newdat <-  gen_new.frame(dat = ic2surv, timepoints = timepoints)

post <- posterior_linpred(m2.stan_gam, newdata = newdat)
pred_surv <- pred_surv(post = post, longdata = newdat)
newdat$surv_mean <- as.vector( apply(pred_surv, 1, mean) ) 

probs <- get_survProb(newdat = newdat)

brier <- pec::pec(probs, Surv(time, status) ~ 1,
                  data = ic2surv,
                  maxtime = max(timepoints),
                  exact = FALSE,
                  exactness = 99L)


#---- Cox model

mod <- coxph(Surv(time, status) ~ chemotherapy*radio_therapy, data = ic2surv)
probs2 <- pec::predictSurvProb(mod,
                              newdata = ic2surv,
                              times = timepoints)
suppressMessages(brier2 <- pec::pec(probs, Surv(time, status) ~ 1,
                                   data = ic2surv,
                                   maxtime = max(timepoints),
                                   exact = FALSE,
                                   exactness = 99L))
