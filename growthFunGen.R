library(scam)
library(tidyverse)

lamean = read.table("./Inputs/LA_MeanSD.txt",header=TRUE)

lamean$sex = as.factor(lamean$sex)

lamean = lamean |>
    filter(AgeClass < 32)

growth_model = scam(mean~s(AgeClass,bs="miso",by=sex),data=lamean,weights=num,gamma=0.3)

lamean$CV = lamean$sd/lamean$mean
sd_model = gam(log(CV)~s(mean,by=sex,k=6),data=lamean,weights=num,gamma=1.5)

library(MASS)
cv_lm = rlm(CV~mean*sex,weights=num,data=lamean,maxit=50)

ggplot(lamean,aes(x=mean,y=CV)) + geom_point(aes(size=sqrt(num),color=sex))

sdpred =expand.grid(mean=0:300,sex=1:2)
sdpred$pred = predict(sd_model,newdata=sdpred)
sdpred$feen = sdpred$mean


sdpredr =expand.grid(mean=0:300,sex=as.factor(1:2))
sdpredr$pred = predict(cv_lm,newdata=sdpredr)
sdpredr$feen = sdpred$mean


ggplot(lamean,aes(x=mean,y=CV)) + geom_point(aes(size=sqrt(num),color=sex)) + geom_line(aes(x=feen,y=pred,color=as.factor(sex)),sdpredr)


pdat1 = expand.grid(AgeClass=seq(0,50,by=0.01),sex=1)
pdat1$pred =  predict(growth_model,newdata=pdat1)
pdat1$epred = pdat$pred

pdat2 = expand.grid(AgeClass=seq(0,50,by=0.01),sex=2)
pdat2$pred =  predict(growth_model,newdata=pdat2)
pdat2$epred = pdat$pred


len_at_age1 = RTMB::splinefun(x=pdat1$AgeClass,y=pdat1$pred)
len_at_age2 = RTMB::splinefun(x=pdat2$AgeClass,y=pdat2$pred)


sdpredr1 =expand.grid(mean=seq(0,300,by=0.01),sex=as.factor(1))
sdpredr1$pred = predict(cv_lm,newdata=sdpredr1)
sdpredr1$feen = sdpredr1$mean
sdpredr1$sd = sdpredr1$pred*sdpredr1$mean

sdpredr2 =expand.grid(mean=seq(0,300,by=0.01),sex=as.factor(2))
sdpredr2$pred = predict(cv_lm,newdata=sdpredr1)
sdpredr2$feen = sdpredr2$mean
sdpredr2$sd = sdpredr2$pred*sdpredr2$mean

###Warning this can go negative......
sd_at_len1 = RTMB::splinefun(x=sdpredr1$mean,y=sdpredr2$sd)
sd_at_len2 = RTMB::splinefun(x=sdpredr1$mean,y=sdpredr2$sd)
