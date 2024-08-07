library(scam)
library(tidyverse)
## generate mean/sd from otolith data
la_dat <- read.table(file = "Inputs/la_data_AC_2010_080224.txt",header = T)

la_summ <- la_dat %>% 
  filter(Sex != 0) %>% 
  group_by(Sex,Age) %>% 
  summarise(mean=mean(Length),
            sd=sd(Length),
            num=n()) %>% 
  filter(!is.na(sd))

ggplot(la_summ,aes(x = Age,y = mean,color=as.factor(Sex)))+geom_point()

## read in mean/sd
#lamean = read.table("./Inputs/LA_MeanSD.txt",header=TRUE)
lamean <- la_summ
lamean$Sex = as.factor(lamean$Sex)

lamean = lamean |>
    filter(Age < 32)

growth_model = scam(mean~s(Age,bs="miso",by=Sex),data=lamean,weights=num,gamma=0.3)

lamean$CV = lamean$sd/lamean$mean
sd_model = gam(log(CV)~s(mean,by=Sex,k=6),data=lamean,weights=num,gamma=1.5)

library(MASS)
cv_lm = rlm(CV~mean*Sex,weights=num,data=lamean,maxit=50)

ggplot(lamean,aes(x=mean,y=CV)) + geom_point(aes(size=sqrt(num),color=Sex))

sdpred =expand.grid(mean=0:300,Sex=1:2)
sdpred$pred = predict(sd_model,newdata=sdpred)
sdpred$feen = sdpred$mean


sdpredr =expand.grid(mean=0:300,Sex=as.factor(1:2))
sdpredr$pred = predict(cv_lm,newdata=sdpredr)
sdpredr$feen = sdpred$mean


ggplot(lamean,aes(x=mean,y=CV)) + geom_point(aes(size=sqrt(num),color=Sex)) + geom_line(aes(x=feen,y=pred,color=as.factor(Sex)),sdpredr)


pdat1 = expand.grid(Age=seq(0,50,by=0.01),Sex=1)
pdat1$pred =  predict(growth_model,newdata=pdat1)
pdat1$epred = pdat1$pred

pdat2 = expand.grid(Age=seq(0,50,by=0.01),Sex=2)
pdat2$pred =  predict(growth_model,newdata=pdat2)
pdat2$epred = pdat2$pred


len_at_age1 = RTMB::splinefun(x=pdat1$Age,y=pdat1$pred)
len_at_age2 = RTMB::splinefun(x=pdat2$Age,y=pdat2$pred)


sdpredr1 =expand.grid(mean=seq(0,300,by=0.01),Sex=as.factor(1))
sdpredr1$pred = predict(cv_lm,newdata=sdpredr1)
sdpredr1$feen = sdpredr1$mean
sdpredr1$sd = sdpredr1$pred*sdpredr1$mean

sdpredr2 =expand.grid(mean=seq(0,300,by=0.01),Sex=as.factor(2))
sdpredr2$pred = predict(cv_lm,newdata=sdpredr1)
sdpredr2$feen = sdpredr2$mean
sdpredr2$sd = sdpredr2$pred*sdpredr2$mean

###Warning this can go negative......
sd_at_len1 = RTMB::splinefun(x=sdpredr1$mean,y=sdpredr2$sd)
sd_at_len2 = RTMB::splinefun(x=sdpredr1$mean,y=sdpredr2$sd)
