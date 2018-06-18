#http://www.rprogramming.info/search/bayesian-lasso/6
# Compare with GP ------------------------------

st_mod2 <- rstanarm::stan_gamm4(status~1+offset(log_dtime)+s(time) + factor(chemotherapy) , data = long_ic2dat , family='poisson')


st_mod <- stan_gamm4(status~1+offset(log_dtime)+s(time, bs='gp', m=c(2))+factor(chemotherapy), data = long_ic2dat , family='poisson')
summary(st_mod)

#Note that the matrix of coefficients includes sigma values, so remove these
preds1 <- as.matrix(st_mod)[,1:13] %*% t(as.matrix(st_mod$x))
preds2 <- as.matrix(st_mod2)[,1:11] %*% t(as.matrix(st_mod2$x))

plot.list2 <- lapply(1:nrow(preds2), function(i){
  long_ic2dat$loghaz <- preds2[i, ]
  dat <- long_ic2dat %>%
    group_by(sample_id) %>%
    mutate(surv = exp( - cumsum(exp(loghaz))) ) %>%
    filter(row_number()==n()) %>%
    ungroup() %>%
    arrange(time) %>%
    select(surv)
}
)
  

plot.matrix <- do.call(cbind, plot.list2)
plot.matrix <- as.matrix(plot.matrix)
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
                          lower = 0, upper = 0)
plot.frame <- rbind(plot.frame, zeros.frame)
mle.surv <- survival::survfit(Surv( time , status) ~ chemotherapy,
                              data = ic2surv  )
obs.mortality <- data.frame(time = mle.surv$time,
                            surv = mle.surv$surv,
                            strata = summary(mle.surv, censored=T)$strata)
zeros <- data.frame(time=0, surv=1, strata=unique(obs.mortality$strata))
obs.mortality <- rbind(obs.mortality, zeros)


ggplot2::ggplot(plot.frame, aes(time, postmean, group = strata)) + 
  geom_line(col = "red", alpha = 0.5) +
  geom_ribbon(aes( ymin = lower,
                   ymax = upper,
                   group = strata), alpha = 0.2, fill = "red") +
  geom_step(data=obs.mortality, aes(time,
                                    surv,
                                    group = strata), col="black", alpha = 0.5)





preds <- list(preds1, preds2) %>%
  map(function(x) {
    x %>%
      t() %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(time = long_ic2dat$time,
             sample_id = long_ic2dat$sample_id) %>%
      gather("sample_id", "loghaz", -time) %>%
      filter(sample_id %in% sample(unique(sample_id), 100))
  }) %>%
  bind_rows(.id="Smooth") %>%
  mutate(Smooth = if_else(Smooth==1, "Gaussian Process (s(time, bs='gp', m=c(2)))", "Thin Plate (s(time))")) %>%
  group_by(sample_id) %>%
  mutate(surv = exp( - cumsum(exp(loghaz))) ) %>%
  filter(row_number()==n()) %>%
  ungroup() %>%
  arrange(time)


mle.surv <- survival::survfit(Surv( time , status) ~ chemotherapy,
                              data = ic2surv  )
obs.mortality <- data.frame(time = mle.surv$time,
                            surv = mle.surv$surv,
                            strata = summary(mle.surv, censored=T)$strata)
zeros <- data.frame(time=0, surv=1, strata=unique(obs.mortality$strata))
obs.mortality <- rbind(obs.mortality, zeros)

ggplot(preds, aes(x=time, y=surv)) +
  geom_step(data=obs.mortality, aes(time,
                                    surv,
                                    group = strata), col="black", alpha = 0.5) +
  geom_line(alpha=0.2, mapping=aes(group=Smooth, col=Smooth)) +
  scale_color_viridis(discrete=TRUE, begin=0, end=0.9) +
  labs(title="Bayesian Modeling of Global Temperature Series",
       subtitle="Comparison of gaussian process and thin-plate splines to\ndeal with temporal autocorrelation using rstanarm",
       caption="100 posterior sample smooths from each model shown\nBased on post by @ucfagls at https://goo.gl/vTRCxB\nData from http://www.metoffice.gov.uk/hadobs/hadcrut4/") +
  ylab("Global Mean Temperature Anomaly") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() + theme(legend.position=c(0.75, 0.75))
