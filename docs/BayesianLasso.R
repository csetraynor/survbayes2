#http://www.rprogramming.info/search/bayesian-lasso/6

library(rstanarm)
library(tidyverse)
library(viridis)

gtemp <- read.delim("http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.5.0.0.annual_ns_avg.txt",
                    sep="", header = FALSE) %>%
  select(1,2) %>%
  set_names(c("Year", "Temperature"))

st_mod2 <- stan_gamm4(Temperature ~ s(Year), data = gtemp)
st_mod <- stan_gamm4(Temperature ~ s(Year, bs='gp', m=c(2)), data = gtemp)

#Note that the matrix of coefficients includes sigma values, so remove these
preds1 <- as.matrix(st_mod)[,1:12] %*% t(as.matrix(st_mod$x))
preds2 <- as.matrix(st_mod2)[,1:10] %*% t(as.matrix(st_mod2$x))

preds <- list(preds1, preds2) %>%
  map(function(x) {
    x %>%
      t() %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(Year = gtemp$Year) %>%
      gather("samp", "Temperature", -Year) %>%
      filter(samp %in% sample(unique(samp), 100))
  }) %>%
  bind_rows(.id="Smooth") %>%
  mutate(Smooth = if_else(Smooth==1, "Gaussian Process (s(Year, bs='gp', m=c(2)))", "Thin Plate (s(Year))"))

ggplot(ppreds, aes(x=Year, y=Temperature)) +
  geom_line(data=gtemp, col="black", alpha = 0.5) +
  geom_line(alpha=0.2, mapping=aes(group=paste(samp, Smooth), col=Smooth)) +
  scale_color_viridis(discrete=TRUE, begin=0, end=0.9) +
  labs(title="Bayesian Modeling of Global Temperature Series",
       subtitle="Comparison of gaussian process and thin-plate splines to\ndeal with temporal autocorrelation using rstanarm",
       caption="100 posterior sample smooths from each model shown\nBased on post by @ucfagls at https://goo.gl/vTRCxB\nData from http://www.metoffice.gov.uk/hadobs/hadcrut4/") +
  ylab("Global Mean Temperature Anomaly") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() + theme(legend.position=c(0.3, 0.75))