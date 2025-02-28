#load libraries
library(offarray)
library(offartmb)
library(tidyverse)

# array spot where Prob = 0
s1 = 2
y1 = 10
lc1 = 9
a1 = 3
b2 = 2

# array spot where Prob = NA
s1 = 2
y1 = 10
lc1 = 2
a1 = 28
b2 = 2

###load TRO_sy, make_fecundity function, and shape/scale tables
# centered around 150cm
make_fecundity <- function(s,len){
  (len/150)^bexp[s]
}

#fecundity functions
la_key <- read.table("Inputs/LA_MeanSD.txt",header = T, sep = "\t",stringsAsFactors = F)

#split la_key into arrays of means and sds
raw_la_means_SA <- la_key %>% 
  arrange(sex,AgeClass) %>% 
  dplyr::select(sex,AgeClass,mean) %>% 
  filter(AgeClass <= 30) %>% 
  spread(AgeClass,mean) %>% 
  dplyr::select(-sex)

raw_la_sd_SA <- la_key %>% 
  arrange(sex,AgeClass) %>% 
  dplyr::select(sex,AgeClass,sd) %>% 
  filter(AgeClass <= 30) %>% 
  spread(AgeClass,sd) %>% 
  dplyr::select(-sex)

## making length/age shape/scale for gamma function
la_means_SA <- offarray(as.matrix(raw_la_means_SA),dimseq=list(s=SEXES, a=A))
la_sd_SA <- offarray(as.matrix(raw_la_sd_SA),dimseq=list(s=SEXES, a=A))

la_shape_SA <- (la_means_SA/la_sd_SA)^2
la_scale_SA <- la_means_SA/la_shape_SA

### stuff that is inside the autoloop #####
    # length of parent
    l1 <- Lvec[lc1]
    # quantile of parent at a1 given l1 length
    qq1 <- pgamma(l1,shape = la_shape_SA[s1,a1],scale = la_scale_SA[s1,a1])
    qq1
    #age of parent when off is born
    a1_at_B2 <- a1 - (y1 - b2)
    
    #length of parent when offspring is born
    l1_at_B2 <- qgamma(qq1,shape = la_shape_SA[s1,a1_at_B2 |> clamp(A)],
                       scale = la_scale_SA[s1,a1_at_B2 |> clamp(A)])
    l1_at_B2
    # generate probability given fecundity of par
    # and TRO at the offspring brith year
    Prob <- (y1 >= b2) * # otherwise Molly was dead before Dolly born
      (a1_at_B2 >= 2) *
      make_fecundity(s1,l1_at_B2) * recip_TRO_SY[s1,b2]
    
    Prob
