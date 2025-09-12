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
Lvec <- seq(10,220,10)
SEXES <- 1:2
LENGTH_CLASSES <- 1:length(Lvec)

## estimate mean/sd for each age class
MeanLength <- lengthage %>% 
  filter(OtoAge <= max(A)) %>% 
  group_by(ObsSex,OtoAge) %>% 
  summarise(Means=mean(Length),SDs=sd(Length)) %>% 
  arrange(OtoAge,ObsSex)

MeanLength$Shape <- (MeanLength$Means/MeanLength$SDs)^2
MeanLength$Scale <- MeanLength$SDs^2/MeanLength$Means


MeanLength_SA <- offarray(as.matrix(MeanLength$Means),dimseq=list(s=SEXES, a=A))
SDLength_SA <- offarray(as.matrix(MeanLength$SDs),dimseq=list(s=SEXES, a=A))
Shape_SA <- offarray(as.matrix(MeanLength$Shape),dimseq=list(s=SEXES, a=A))
Scale_SA <- offarray(as.matrix(MeanLength$Scale),dimseq=list(s=SEXES, a=A))

## perform back-calculations for all the quantile estimates
## make offarray the size of s1,a1,lc1
LenatB2_SYLAB <- autoloop( 
  s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, a1=A,
  b2=POPY, {
    # length of parent
    l1 <- Lvec[lc1]
    # quantile of parent at a1 given l1 length
    qq1 <- pgamma(l1,shape = Shape_SA[s1,a1],scale = Scale_SA[s1,a1])+1e-15
    qq1 <- (a1 < 5) * pgamma(l1,shape = Shape_SA[s1,a1],scale = Scale_SA[s1,a1])-1e-15
    qq1 <- (a1 >= 20) * pgamma(l1,shape = Shape_SA[s1,a1],scale = Scale_SA[s1,a1])+1e-15
    #age of parent when off is born
    a1_at_B2 <- a1 - (y1 - b2)
    
    #length of parent when offspring is born
    l1_at_B2 <- qgamma(qq1,shape = Shape_SA[s1,a1_at_B2 |> clamp(A)],
                       scale = Scale_SA[s1,a1_at_B2 |> clamp(A)])
    (l1_at_B2 <= l1) *
    l1_at_B2
  })
