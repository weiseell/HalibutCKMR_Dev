## creating length quantile bins based on otolith age data

## load data
lengthage <- read.csv("E:/ClonedRepositories/ForTheHalibut/Input/Armsworthy&Campana_ladata_nofilter.csv")

## load libraries
library(tidyverse)
library(offarray)
library(mvbutils)
## load current model constants
A <- 2:30
POPY <- 1:15
SAMPY <- 10:15
Lvec <- 1:10
SEXES <- 1:2

## estimate mean/sd for each age class
MeanLength <- lengthage %>% 
  filter(OtoAge <= max(A)) %>% 
  group_by(ObsSex,OtoAge) %>% 
  summarise(Means=mean(Length),SDs=sd(Length))

## calculate representative quantiles for each age class
LengthQ <- data.frame(matrix(data = NA,nrow = nrow(MeanLength),ncol = max(Lvec)))
for (i in 1:nrow(MeanLength)) {
  LengthQ[i,] <- qnorm(seq(0.05,0.95,by=0.1),mean = MeanLength$Means[i],sd = MeanLength$SDs[i])
}

LengthQ <- cbind(MeanLength,LengthQ)
LengthQ <- LengthQ %>% dplyr::select(-Means,-SDs)

# turn lengthQ into an array
LengthQ <- LengthQ %>% 
  gather(key = "lc", value = "lengthQ",-ObsSex,-OtoAge) %>% 
  separate(lc,into = c("tmp","lc"),sep = "X") %>% 
  dplyr::select(-tmp) %>% 
  mutate(lc=as.numeric(lc)) %>% 
  arrange(ObsSex,OtoAge,lc)

LengthQ_SAL <- offarray(x=0,dimseq= list(s1=SEXES,a1=A,lc1=Lvec))
for (i in 1:nrow(LengthQ)) {
  LengthQ_SAL[LengthQ$ObsSex[i],LengthQ$OtoAge[i],LengthQ$lc[i]] <- LengthQ$lengthQ[i]
}

## perform back-calculations for all the quantile estimates
## make offarray the size of s1,a1,lc1
LengthQ_SYLAB <- autoloop( 
  s1=SEXES, y1=SAMPY, lc1=Lvec, a1=A,
  b2=POPY, {
    # length of parent at catch
    l1 <- LengthQ_SAL[SLICE=s1,SLICE=a1,SLICE=lc1]
    
    # backcalculate parent age at birth year
    a1_at_B2 <- a1 - (y1 - b2)
    l1_at_B2 <- LengthQ_SAL[s1,a1_at_B2 |> clamp(A),lc1]
    l1_at_B2
  })
