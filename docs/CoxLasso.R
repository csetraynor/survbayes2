#mc-cross-validation.R

#Load libraries
library(dplyr)

#Load data
data("ic2surv")
data("iclust2_glmnet")

#Extract features, see my_replace in utils
iclust2_features <- extract_features(iclust2_glmnet)
iclust2_features$feature <- iclust2prog::my_replace(iclust2_features$feature)



#Get brca iclust2

iclust2 <- ic2surv[, c(iclust2_features$feature, "radio_therapy", "chemotherapy", "time", "status")]


model_fit <- function(dat){

# create predictor matrix
x <- dat %>% dplyr::select(-time, -status) %>%
  dplyr::mutate(chemotherapy = chemotherapy == "YES",
                radio_therapy = radio_therapy == "YES")

###### Fit models
##################################################
#grouped enet
p.fac = rep(1, ncol(x))
p.fac[match(c("npi","age_std", "radio_therapy", "chemotherapy"), colnames(x))] = 0
#prepare
x <- as.matrix(x)
y <- as.matrix(dat %>%
                 dplyr::select(time = time,
                               status = status), ncol = 2)

### Apply elastic net with a=1 to perform Lasso
mod <-  glmnet::cv.glmnet(x, y, family = "cox",
                          grouped = TRUE,
                          alpha = 1,
                          parallel = FALSE, penalty.factor = p.fac)
return(mod)
}

mod <- model_fit(iclust2)
ic2_lasso_gwt <- mod
devtools::use_data(ic2_lasso_gwt)
saveRDS(mod, "iclust2_lasso.RDS")
rm(mod)
rm(iclust2)
modE <- model_fit(ERpos)
saveRDS(modE, "erpos_lasso.RDS")


