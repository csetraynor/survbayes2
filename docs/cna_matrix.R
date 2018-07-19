options(expressions = 5e5)
library(dplyr)

t_cna_matrix <- readRDS("C:/RFactory/parallel/t_cna_matrix.RDS")

t_cna_matrix <- readRDS("/home/mtr/rfactory/predsurv/vignettes/t_cna_matrix.RDS")

t_cna_matrix <- as.data.frame(t_cna_matrix)

t_cna_matrix <- t_cna_matrix %>% mutate_all(as.factor)

X <- model.matrix( ~ ., t_cna_matrix[,1:10])[,-1]

X <- as.data.frame(X)

for(i in 11:ncol(t_cna_matrix)){
  
  #In case of more than one level
  if(nlevels(t_cna_matrix[,i]) > 1){
    coluna <- as.data.frame(t_cna_matrix[,i])
    colnames(coluna) <- colnames(t_cna_matrix[i])
    matrixna <- model.matrix(~ ., coluna)[,-1]
    matrixna <- as.data.frame(matrixna)
    print(colnames(matrixna))
    X <- cbind(X, matrixna)
  }
  #Transform to numeric
  else{
    coluna <- as.numeric(as.factor(t_cna_matrix[,i]))
    coluna <- as.data.frame(coluna)
    colnames(coluna) <- colnames(t_cna_matrix[i])
    X <-cbind(X, coluna)
  }
}


##Check zero var
t_cna_expression <- as.data.frame(X)
test <- zeroVar(t_cna_expression)
if(length(test) > 0){
  t_cna_expression <- t_cna_expression[,-test]
}

saveRDS(t_cna_expression, "/home/mtr/rfactory/cna_expression.RDS")


Y <- readRDS("/home/mtr/rfactory/brca_data.RDS")
### cbind clinical and genomic data
brca_data <- cbind(Y, t_cna_expression)

Y <- Y[,-grep("cna", colnames(brca_data))]

#Load libraries
library(dplyr)
require(doMC)

#Get brca iclust2
iclust2 <- brca_data[brca_data$intclust == 2,]
iclust2$intclust <- NULL
brca_data$intclust <- NULL
hold_out <- lapply(mc_samples, function(x)
  iclust2$patient_id[x])




#Train : Elastic net model
rocky <- function(holdout){
  #### obtain training datset
  if('patient_id' %in% colnames(survdata)){
    train <- survdata[!( (survdata %>%
                            dplyr::select(patient_id) %>%
                            unlist) %in%
                           (holdout) ),];
    #unselect subject
    train <- train %>% dplyr::select(- patient_id)
  }
  
  # create predictor matrix
  x <- train %>% dplyr::select(-os_months, -os_deceased)
  
  ###### create fold id for CV in glmnet
  set.seed(9)
  foldid <-  caret::createFolds(train %>% select(os_deceased) %>% unlist,
                                k = 10, list = FALSE)
  ###### Fit models
  ##################################################
  #grouped enet
  p.fac = rep(1, ncol(x))
  p.fac[match(c("npi","age_std"), colnames(x))] = 0.5
  #prepare
  x <- as.matrix(x)
  y <- as.matrix(train %>%
                   dplyr::select(time = os_months,
                                 status = os_deceased), ncol = 2)
  #require(doMC)
  #registerDoMC(cores=nw)
  ### Apply elastic net with a=0.8 closer to Lasso
  mod <-  glmnet::cv.glmnet(x, y, family = "cox",
                            grouped = TRUE,
                            lambda.min.ratio = lambda, alpha = 0.8,
                            foldid = foldid,
                            parallel = FALSE, penalty.factor = p.fac)
  # find optimised lambda
  optimal <- as.matrix(coef(mod, s = "lambda.min"))
  optimal <- as.data.frame(optimal)
  colnames(optimal) <- "mod"
  optimal$Hugo_symbol <- rownames(optimal)
  optimal <-  optimal %>% filter(mod != 0)
  return(optimal)
}

# This assumes that length(datlist) is much less than ncores
ncores <-  parallel::detectCores()
m <- 100
#nw <- ncores %/% m

survdata = iclust2
lambda = 0.0001
system.time(iclust2_fit <- mclapply(hold_out, rocky, mc.cores=m))


saveRDS(iclust2_fit, "cv_iclust2_fit_lasso_cna_10.RDS")

survdata = brca_data
lambda = 0.001
system.time(pooled_fit <- mclapply(hold_out, rocky, mc.cores=m ) )

saveRDS(pooled_fit, "cv_pooled_fit_lasso_cna_10.RDS")

#10
