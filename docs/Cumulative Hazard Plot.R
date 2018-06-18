#---------Cumulative Hazard Plot----------
library(plyr)
library(survival)
library(flexsurv)
library(ggplot2)    # for ggplot function
library(gridExtra)  # for grid arrange function


#-------Load NONMEM file


## Name of working directory
wdir <- "C:/Users/carlos.traynor/Documents/Project/BPHrepsnew/NONMEM TTE files"

fileName <- "sdtab003"

d <- read.table(paste(wdir,"/",fileName,sep=""),header=T,skip=1)

#---Load Observed Data

obs <- read.csv('C:/Users/carlos.traynor/Documents/Project/BPHrepsnew/Dataset/_first_treatment_time_to_first_event_dataset.csv', header = T, sep = ",")


#-----Load Covariates

fileCov <- "cotab003"
c <-  read.table(paste(wdir,"/",fileCov,sep=""),header=T,skip=1)

p <- join(c,d)



#----Preparation Observed Data

obs <- obs[obs$MDV==0,]
obs$TREAT <- c("Placebo","Tamsulosin","Dutasteride","Watchful Waiting","Combination")[obs$TREAT]
obs$SMOKE <- c("NO","YES","NA")[obs$SMOKE]
obs$DV <- obs$DV + 1

#------------------------------Preparation NONMEM

p$SMOKE <- c("NO","YES","NA")[p$SMOKE]

#----Confidence of intervals---#
#---Load Covariance Matrix----#
var <- read.csv("C:/Users/carlos.traynor/Documents/Project/BPHrepsnew/NONMEM TTE files/covariance matrix model3.csv")
v <- rep(0, dim(var)[1])

j <- 2
for (ii in 2:dim(var)[1]){
  v[j-1] <-var[j-1,j]
  j <- j+1
}

se <- sqrt(v)

write.csv(se, "C:/Users/carlos.traynor/Documents/Project/BPHrepsnew/Dataset/se.csv")

p$A1 <-0
p$A1[p$BIPSSGRP2==3] <- 1

p$A2 <- 0
p$A2[p$BPSAGRP==2] <- 1

p$A3 <- 0
p$A3[p$BPSAGRP==3] <- 1

p$A4 <- 0
p$A4[p$SMOKE=="YES"] <- 1

p$T1 <- 0
p$T1[p$TREAT==2] <- 1

p$T2 <- 0
p$T2[p$TREAT==3] <- 1 

p$T3 <- 0
p$T3[p$TREAT==4] <- 1

p$T4 <- 0
p$T4[p$TREAT==5] <- 1

parameters <- read.csv("C:/Users/carlos.traynor/Documents/Project/BPHrepsnew/NONMEM TTE files/parameters run 3.csv", 
                       sep=",", row.names=NULL)
par <- parameters[3:11,2:3]

#95%confidence interval change to 1.96 for 95%
v <- as.numeric(as.character(par[,2]))
upar <- v[1:8]+1.645*se[1:8]
lpar <- v[1:8]-1.645*se[1:8]


p$URS <- upar[1]*p$A1+upar[2]*p$A2+upar[3]*p$A3+upar[4]*p$A4+upar[5]*p$T1+upar[6]*p$T2+upar[7]*p$T3+upar[8]*p$T4
p$LRS<- lpar[1]*p$A1+lpar[2]*p$A2+lpar[3]*p$A3+lpar[4]*p$A4+lpar[5]*p$T1+lpar[6]*p$T2+lpar[7]*p$T3+lpar[8]*p$T4
p$RS <- v[1]*p$A1+v[2]*p$A2+v[3]*p$A3+v[4]*p$A4+v[5]*p$T1+v[6]*p$T2+v[7]*p$T3+v[8]*p$T4

ulambda <- (v[9]+1.96*se[9])/10000
llambda <- (v[9]-1.96*se[9])/10000
p$HAZ0 <- ((v[9]*exp(p$RS))/10000)*p$TIME
p$UL <- ulambda*exp(p$URS)*p$TIME
p$LL <- llambda*exp(p$LRS)*p$TIME



base <- p[p$BPSAGRP==3& ( p$BIPSSGRP2 ==1 ) & p$SMOKE=="NO" &p$MDV==0,]
base$STRATA <- c("Placebo","Tamsulosin","Dutasteride","Watchful Waiting","Combination")[base$TREAT]

base <- base[c("TIME","SUR","HAZ0","UL","LL","STRATA")]
#---Observeved ----->
base.fit <- survfit(Surv(TIME, DV)~TREAT  , data=obs[obs$BPSAGRP==3 & (obs$BIPSSGRP2 ==1 | obs$BIPSSGRP2 ==2) & obs$SMOKE=="NO",])
obs.mortality <- data.frame(TIME=base.fit$time, SUR=base.fit$surv, STRATA=summary(base.fit, censored=T)$strata, UL=base.fit$surv+1.96*base.fit$std.err, LL=base.fit$surv-1.96*base.fit$std.err)
zeros <- data.frame(TIME=0, SUR=1, STRATA=unique(obs.mortality$STRATA), UL=1 , LL=1)
obs.mortality <- rbind(obs.mortality, zeros)
obs.mortality$HAZ0 <- -log(obs.mortality$SUR) # plot survival function or hazard function
obs.mortality$UL <- -log(obs.mortality$UL)
obs.mortality$LL <- -log(obs.mortality$LL)
obs.mortality$STRATA <- as.factor(obs.mortality$STRATA)
levels(obs.mortality$STRATA) <- c("Combination","Dutasteride","Placebo","Tamsulosin","Watchful Waiting")



#----join---
new <- rbind(transform(base, SOURCE="Prediction"),transform(obs.mortality, SOURCE="Observation"))

#---Plot GOF----#
ggplot(new[new$TIME<1460,], aes(TIME, HAZ0,  linetype = factor(SOURCE)))+geom_step()+
  labs(title= "Goodness Of Fit", x="time[days]", y="Cumulative Hazard")+facet_wrap(~STRATA)+geom_ribbon(aes(ymin=LL, ymax=UL), alpha=1/3)




#--------------------------------------------------------------End OF Script


# 
# 
# 
# ggplot(base[base$TIME<1460,], aes(TIME, H0))+geom_step(aes(col="VPC"))+
#          geom_ribbon(aes(ymin=LH0, ymax=UH0), alpha=1/6)+
#          geom_step(data=new[new$TIME<1460&new$SOURCE=="observation",], aes(TIME, HAZ0, col="Orig"), linetype=2)+labs(col="")+
#          facet_wrap(~STRATA)+labs(title="Model Prediction for Baseline Characteristics", x="TIME[days]", y="Cumulative Hazard")
# 
# 
# #---For PSAGRP3
# 
# psa <- p[p$BPSAGRP==3 & ( p$BIPSSGRP2 ==1 | p$BIPSSGRP2 ==2) & p$SMOKE=="NO",]
# 
# psa$STRATA <- c("Placebo","Tamsulosin","Dutasteride","Watchful Waiting","Combination")[psa$TREAT]
# 
# ggplot(psa[psa$TIME<1460,], aes(TIME, H0))+geom_step(aes(col="VPC"))+
#   geom_ribbon(aes(ymin=LH0, ymax=UH0), alpha=1/6)+
#   geom_step(data=new[new$TIME<1460&new$SOURCE=="observation",], aes(TIME, HAZ0, col="Orig"), linetype=2)+labs(col="")+
#   facet_wrap(~STRATA)+labs(title="IC 90%", x="TIME[days]", y="Cumulative Hazard")
# 
# bph$BPSA <- as.numeric(as.character(bph$BPSA))
# bph$BPSA[bph$BPSA =='NaN'] <- NA
# bph$BPSA[bph$BPSA ==''] <- NA
# 
# #------Plot Survival Curves for Baseline Covariates 
# 
# 
# 
# #---Observeved ----->
# base.fit <- survfit(Surv(TIME, DV)~TREAT  , data=obs[obs$BPSAGRP==3 & (obs$BIPSSGRP2 ==1 | obs$BIPSSGRP2 ==2) & obs$SMOKE=="NO",])
# obs.mortality <- data.frame(TIME=base.fit$time, SUR=base.fit$surv, STRATA=summary(base.fit, censored=T)$strata, UL=base.fit$surv+1.96*bph.fit$std.err, LL=base.fit$surv-1.96*bph.fit$std.err))
# zeros <- data.frame(TIME=0, SUR=1, STRATA=unique(obs.mortality$STRATA), UL=1 , LL=1)
# obs.mortality <- rbind(obs.mortality, zeros)
# obs.mortality$HAZ0 <- -log(obs.mortality$SUR) # plot survival function or hazard function
# 
# 
# #---Predicted------>
# d2 <- p[p$BPSAGRP==1 & ( p$BIPSSGRP2 ==1 | p$BIPSSGRP2 ==2) & p$SMOKE=="NO",c('TIME', 'SUR','TREAT')]
# d2$HAZ0 <- -log(d2$SUR)
# d2$STRATA <- as.factor(as.character(d2$TREAT))
# d2$STRATA <- c("TREAT=Placebo","TREAT=Tamsulosin","TREAT=Dutasteride","TREAT=Watchful Waiting","TREAT=Combination")[d2$STRATA]
# d2$TREAT <- NULL
# 
# library(OIsurv) 
# confBands(base.fit, confLevel=0.95, type="hall")
# 
# 
# #----join---
# new <- rbind(transform(d2, SOURCE="prediction"),transform(obs.mortality, SOURCE="observation"))
# new$STRATA <- as.factor(new$STRATA)
# levels(new$STRATA) <- c("Combination","Dutasteride","Placebo","Tamsulosin","Watchful Waiting")
# 
# #---Plot GOF----#
# ggplot(new[new$TIME<1460,], aes(TIME, HAZ0, colour=factor(STRATA), linetype = factor(SOURCE)))+geom_step()+
#   labs(title= "GOF Baseline Explanatory Variables", x="time[days]", y="Cumulative incidence[%]")+facet_wrap(~STRATA)
# 
# 
# 
# #---Observeved ----->
# base.fit <- survfit(Surv(TIME, DV)~TREAT  , data=obs[obs$BPSAGRP==3 & (obs$BIPSSGRP2 ==1|obs$BIPSSGRP2==2 ) & obs$SMOKE=="NO",])
# obs.mortality <- data.frame(TIME=base.fit$time, SUR=base.fit$surv, STRATA=summary(base.fit, censored=T)$strata)
# zeros <- data.frame(TIME=0, SUR=1, STRATA=unique(obs.mortality$STRATA))
# obs.mortality <- rbind(obs.mortality, zeros)
# obs.mortality$HAZ0 <- -log(obs.mortality$SUR) # plot survival function or hazard function
# 
# 
# #---Predicted------>
# d2 <- p[p$BPSAGRP==3 & ( p$BIPSSGRP2 ==1 | p$BIPSSGRP2 ==2) & p$SMOKE=="NO",c('TIME', 'SUR','TREAT')]
# d2$HAZ0 <- -log(d2$SUR)
# d2$STRATA <- as.factor(as.character(d2$TREAT))
# d2$STRATA <- c("TREAT=Placebo","TREAT=Tamsulosin","TREAT=Dutasteride","TREAT=Watchful Waiting","TREAT=Combination")[d2$STRATA]
# d2$TREAT <- NULL
# 
# library(OIsurv) 
# confBands(base.fit, confLevel=0.95, type="hall")
# 
# 
# #----join---
# new <- rbind(transform(d2, SOURCE="prediction"),transform(obs.mortality, SOURCE="observation"))
# new$STRATA <- as.factor(new$STRATA)
# levels(new$STRATA) <- c("Combination","Dutasteride","Placebo","Tamsulosin","Watchful Waiting")
# 
# #---Plot GOF----#
# ggplot(new[new$TIME<1460,], aes(TIME, HAZ0, colour=factor(STRATA), linetype = factor(SOURCE)))+geom_step()+
#   labs(title= "Goodness Of Fit", x="time[days]", y="Cumulative Hazard")+facet_wrap(~STRATA)

#------------load data----------------
# bph <- read.csv('C:/Users/carlos.traynor/Documents/Project/BPHrepsnew/Dataset/_first_treatment_time_to_first_event_dataset.csv', header = T, na='.')
# #------------------------------------
# 
# #-----Preparation-----#
# bph <- bph[bph$MDV==0,]
# 
# bph$TREAT <- c("1Placebo","2Tamsulosin","Dutasteride","Watchful Waiting","Combination")[bph$TREAT]
# bph$SMOKE <- c("NO","YES","NA")[bph$SMOKE]
# 
# 
# 
# ph.fit <- coxph(Surv(TIME, DV)~as.factor(BIPSSGRP2)+as.factor(BPSAGRP)+as.factor(SMOKE)+as.factor(TREAT), data=bph[bph$TIME>0 &bph$SMOKE!="NA" & bph$BIPSSGRP2 != 4 & bph$BPSAGRP != 4,])
# 
# 
# exp.fit <- survreg(Surv(TIME, DV)~as.factor(BIPSSGRP2)+as.factor(BPSAGRP)+as.factor(SMOKE)+as.factor(TREAT), data=bph[bph$TIME>0 &bph$SMOKE!="NA" & bph$BIPSSGRP2 != 4 & bph$BPSAGRP != 4,], dist = "exp")
# 
# exp.flex <- flexsurvreg(Surv(TIME, DV)~as.factor(BIPSSGRP2)+as.factor(BPSAGRP)+as.factor(SMOKE)+as.factor(TREAT), data=bph[bph$TIME>0 &bph$SMOKE!="NA" & bph$BIPSSGRP2 != 4 & bph$BPSAGRP != 4,], dist = "exp")
# 
# 
# #--------BIPSSGRP 2 contains 0.
# 
# bph$BIPSSGRP2[bph$BIPSSGRP2==1] <- "1-BASE"
# bph$BIPSSGRP2[bph$BIPSSGRP2==2] <- "1-BASE"
# 
# 
# 
# ph.fit <- coxph(Surv(TIME, DV)~as.factor(BIPSSGRP2)+as.factor(BPSAGRP)+as.factor(SMOKE)+as.factor(TREAT), data=bph[bph$TIME>0 &bph$SMOKE!="NA" & bph$BIPSSGRP2 != 4 & bph$BPSAGRP != 4,])
# 
# 
# exp.fit <- survreg(Surv(TIME, DV)~as.factor(BIPSSGRP2)+as.factor(BPSAGRP)+as.factor(SMOKE)+as.factor(TREAT), data=bph[bph$TIME>0 &bph$SMOKE!="NA" & bph$BIPSSGRP2 != 4 & bph$BPSAGRP != 4,], dist = "exp")
# 
# exp.flex <- flexsurvreg(Surv(TIME, DV)~as.factor(BIPSSGRP2)+as.factor(BPSAGRP)+as.factor(SMOKE)+as.factor(TREAT), data=bph[bph$TIME>0 &bph$SMOKE!="NA" & bph$BIPSSGRP2 != 4 & bph$BPSAGRP != 4,], dist = "exp")
