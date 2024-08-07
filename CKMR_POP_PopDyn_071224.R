### Just simulated Data input pop dyn test
# random-walk abundance
library(tidyverse)
library(RTMB)
library(offarray)
library(mvbutils)
library(offartmb)

#load simulated data inputs
load("Inputs/CensusSize_Sim_071224.rda")
load("Inputs/Ncomp_POPs_SYLSYL_070924.rda")
load("Inputs/NPOPs_SYLSYL_070924.rda")

#source functions for RTMB model
source("ModelFunction/prob_la.R")

Narray = array(CS_181_195$abundance,c(50,15,2))
dat <- list()
dat$Narray <- Narray
dat$A <- 2:30
dat$POPY <- 1:15
dat$SAMPY <- 10:15
dat$Lvec <- seq(10,220,10)
dat$LENGTH_CLASSES <- 1:length(dat$Lvec)
dat$SEXES <- 1:2
dat$numq <- 5
dat$la_key <- read.table("Inputs/LA_MeanSD.txt",header = T, sep = "\t",stringsAsFactors = F)

#number of POPs and comps as array (from script 2 that processes sim results)
dat$n_POP_SYLSYL <- n_POP_SYLSYL
dat$n_comps_POP_SYLSYL <- n_comps_POP_SYLSYL

#split la_key into arrays of means and sds
raw_la_means_SA <- dat$la_key %>% 
  arrange(sex,AgeClass) %>% 
  select(sex,AgeClass,mean) %>% 
  filter(AgeClass <= 30) %>% 
  spread(AgeClass,mean) %>% 
  select(-sex)

raw_la_sd_SA <- dat$la_key %>% 
  arrange(sex,AgeClass) %>% 
  select(sex,AgeClass,sd) %>% 
  filter(AgeClass <= 30) %>% 
  spread(AgeClass,sd) %>% 
  select(-sex)

dat$la_means_SA <- offarray(as.matrix(raw_la_means_SA),dimseq=list( s=dat$SEXES, a=dat$A))
dat$la_sd_SA <- offarray(as.matrix(raw_la_sd_SA),dimseq=list( s=dat$SEXES, a=dat$A))

# create prob len at age array
raw_prob_len_at_age <- prob_la(Lmax = max(dat$Lvec),
                           ages = dat$A,
                           binsize = 10,
                           dat = dat$la_key)

# Make it into offarray
dat$prob_len_at_age <- offarray(raw_prob_len_at_age, 
        dimseq=list( s=dat$SEXES, lc=dat$LENGTH_CLASSES, a=dat$A))

## create parameter list
parm <- list()
#log recruitment - random effect length popdynyears
parm$log_rec <- log(Narray[1,,1] + Narray[1,,2])
parm$log_rec_sd = log(0.3)
#parm$log_rec <- rep_len(x = 0.1, length.out = dat$Y)
#total log abundance - split into Y/A abundance within the model
#parm$log_init_abundance <- 0.1
parm$log_init_abundance <- log(sum(Narray[-1,1,]))
# average log survival for M and F
parm$log_Z <- log(c(0.2,0.2))
#a and b values for the fecundity equations
#!# fix these for now
parm$log_bexp <- c(0.1,0.1)

new_f <- function(parm) reclasso( by=parm, {
  getAll(dat,parm)
  
  make_fecundity <- function(s,len){
    (len/150)^bexp[s]
  }

  #exp the sd
  rec_sd = exp(log_rec_sd)
  ##An array for male/female fecundity
  bexp <- exp(log_bexp)
  rec <- exp(log_rec)
  ##Array for average fecundity by sex and age
  ## create prob_by_length array and then fecun array
  
  #prob_len_at_age <- autoloop(
  #  s=SEXES, lc=LENGTH_CLASSES, a=AGES,
  #  pnorm( ) - pnorm( )
  #)
  
  fec_sa <- autoloop(
      s=SEXES, a=A, 
      SUMOVER=list(lc=LENGTH_CLASSES), {
    make_fecundity(s,Lvec[lc]) * prob_len_at_age[s,lc,a]
        #((Lvec[lc]/150)^bexp[s]) * prob_len_at_age[s,lc,a]
  })
  
  fec_sa_quant <- autoloop(
    s=SEXES, a=A, qq=(1:nquant), {
      #calculating length for quantile at mean/sd
      l_q <- qnorm(qq/(nquant + 1),mean = la_means_SA[s,a],sd = la_sd_SA[s,a])
      
      #make fecundity for length given prob
      make_fecundity(s,l_q) * (1/nquant)
    })
  
  ##Plus group for later
  N = offarray(0,dimseq = list(SEXES=SEXES,POPY=POPY,AGE=A))

  nll = 0
  ##Random walk for recruitment
  for(y in 2:length(POPY)){
    nll = nll - dnorm(log_rec[y],log_rec[y-1],rec_sd,TRUE)
  }
  
  ##Put in recruitment
  for(s in 1:2){
    N[s,,2] = rec/2
  }
  
  Z = 0 * N
  ##Male Mortality
  Z[1,,] = exp(log_Z[1])
  ##Female
  Z[2,,] = exp(log_Z[2])
  
  cumsurv_m = rev(cumsum(Z[SLICE=1,SLICE=1,head(A,-1)])) 
  cumsurv_f = rev(cumsum(Z[SLICE=2,SLICE=1,head(A,-1)]))
  total_surv = sum(cumsurv_m+cumsurv_f)
  
  tmp <- (cumsurv_f/total_surv)*exp(log_init_abundance)
  tmp <- c( N[SLICE=2,SLICE=1,SLICE=2], tmp)
  N[2,1,] <- tmp 
  
  tmp <- (cumsurv_m/total_surv)*exp(log_init_abundance)
  tmp <- c( N[SLICE=1,SLICE=1,SLICE=2], tmp)
  N[1,1,] <- tmp
  
  #N[1,1,A[-1]] <- (cumsurv_m/total_surv)*exp(log_init_abundance)
  #N[2,1,A[-1]] <- (cumsurv_f/total_surv)*exp(log_init_abundance)
  
  ## add the plus group to year 1
  Abar_plus <- offarray(0,dimseq = list(SEXES=SEXES,POPY=POPY))
  
  Abar_plus[,1] <- max(A) + 1/(1-exp(-Z[,SLICE=1,SLICE=max(A)]))
  
  adults <- A[-1]
  for (y in 2:length(POPY)) {
    ## run abundance for the next year
    N[,y,adults] <- autoloop(s = SEXES, a = adults,{
      N[s,y-1,a-1]*exp(-Z[s,y-1,a-1])
    })
    
    #add the plus group
    N[,y,max(A)] <- N[,SLICE=y,SLICE=max(A)] + N[,SLICE=y-1,SLICE=max(A)] * exp(-Z[,SLICE=y-1,SLICE=max(A)])
    
    #add the Abar equation for the plus group
    Abar_plus[,y] <- c(((Abar_plus[,SLICE=y-1]+1) * N[,SLICE=y-1,SLICE=max(A)] + 
                          max(A) * N[,SLICE=y-1,SLICE=max(A)-1])/
      (N[,SLICE=y-1,SLICE=max(A)] + N[,SLICE=y-1,SLICE=max(A)-1]))
    
  }
  
  ## calculate survival array to deal with plus group
  # used in half-sib probability
  maxage <- max(A) + max(POPY)
  superAge <- 2:maxage
  
  ## run abundance for the next year
  superN <- offarray(0,dimseq = list(SEXES=SEXES,POPY=POPY,AGE=superAge))
  superN[,,A] <- N
  
  superN[,2:max(POPY),max(A)] <- 0
  for (y in 2:length(POPY)) {
    superN[,y,max(A):max(maxage)] <- superN[,y-1,(max(A):max(maxage))-1] * exp(-Z[,y-1,max(A)])
  }
  
  Pr_Surv_SYAY <- autoloop(
    spp=SEXES,y1=POPY,app=A,y2=POPY,{
      superN[spp,y2,app + (y2-y1)]/superN[spp,y1,app]
    }
  )
  
  TRO_SY <- autoloop(
      s=SEXES, y=POPY,
      SUMOVER=list(a=A),{
    N[s,y,a] * fec_sa[s,a]    
  })
  
  recip_TRO_SY <- 1/TRO_SY
  
  ## nll using CKMR data
  ## nll using CKMR data
  #Pr_MOP_SYLASYA <- array(0,c(2,Y,length(Lvec),A,2,Y,A))
  # generating the POP probabilities
  # This does all ages for IDEAl measurement
  # autoloop() both CREATES the array and FILLS IT IN
  # looping over all the "equals" at the start
  # This does the plus-group wrong. We will fix immediately afterwards...
  Pr_POP_SYLAB <- autoloop( 
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, a1=A,
    b2=POPY, {
      #sd that parent length is currently from the mean
      
      l1 <- Lvec[lc1]
      sd1 <- (l1 - la_means_SA[s1,a1])/la_sd_SA[s1,a1]
      
      #age of parent when off is born
      a1_at_B2 <- a1 - (y1 - b2)
      #length of parent when offspring is born
      
      l1_at_B2 <- la_means_SA[s1,a1_at_B2 |> clamp( A)] + 
        sd1 * la_sd_SA[s1,a1_at_B2 |> clamp( A)]
      l1_at_B2 <- ifelse( l1_at_B2<0, 0, l1_at_B2) # avoid min, max, etc; diffable (?)
      # ... preferable to do ifelse() *before* make_fecundity(), so we avoid neg nums to pwr...
      
      #!# switch fecundity input from age-based to length-based in the fecundity
      Prob <- 
        (y1 >= b2) * # otherwise Molly was dead before Dolly born
        (a1_at_B2 >= 2) *
        make_fecundity(s1,l1_at_B2) * recip_TRO_SY[s1,b2]
      
      Prob
    })
  
                # PLUS GROUP FIXUP HERE
                # adjusting probabilities to account for plus group
                # age of parent is 30+ instead of only 30
                
                #Pr_POP_SYLAB_plus <- autoloop(s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES,
                #                              b2=POPY, SUMOVER=list(q=1:numq), {
                                                #sd that parent length is currently from the mean
                  
                  
                  
                  # calculate what a1 is based on the length
                  #a1 <- qexp(q/(numq + 1)) * (Abar_plus[s1,y1] - max(A)) + max(A)
                  
                  #l1 <- Lvec[lc1]
                  #!# CAN FIX THIS IN THE FUTURE BY CHANGING THE MEANS/SD FOR LENGTH/AGE
                  # TO A FUNCTION BASED ON LA DATA
                  
                  #sd1 <- (l1 - la_means_SA[s1,a1])/la_sd_SA[s1,a1]
                  
                  #age of parent when off is born
                  #a1_at_B2 <- a1 - (y1 - b2)
                  #length of parent when offspring is born
                  
                  #l1_at_B2 <- la_means_SA[s1,a1_at_B2 |> clamp( A)] + 
                  #    sd1 * la_sd_SA[s1,a1_at_B2 |> clamp( A)]
                  #Prob <- ifelse(l1_at_B2 > 0,
                  #!# switch fecundity input from age-based to length-based in the fecundity
                  #  (y1 >= b2) * # otherwise Molly was dead before Dolly born
                  #  (a1_at_B2 >= 2) *
                  #  make_fecundity(s1,l1_at_B2) * inv_TRO_SY[s1,b2],
                  #  0)
                  
                  #Prob
                #})
                
    #Pr_POP_SYLAB[,,,,max(A)] <- Pr_POP_SYLAB_plus
                            
    num_Pr_A_SYL <- autoloop(
      a=A, s=SEXES, y=SAMPY, lc=LENGTH_CLASSES,
      prob_len_at_age[s,lc,a] * N[s,y,a]
    )              

    #denom_SYL <- sumover(num_Pr_A_SYL,'a')
    recip_denom_SYL <- 1 / autoloop( 
      indices=dimseq( num_Pr_A_SYL)[-1],   # see num_Pr_A_SYL above for actual list
      SUMOVER=dimseq( num_Pr_A_SYL)[1],    # just the "a"
      num_Pr_A_SYL[a, s, y, lc]
    )
    
    Pr_A_SYL <- autoloop(
      a=A, s=SEXES, y=SAMPY, lc=LENGTH_CLASSES,
      num_Pr_A_SYL[a, s, y, lc] * recip_denom_SYL[s, y, lc]
    )
    
    
    Pr_POP_SYLB <- autoloop(s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES,b2=POPY,
                            SUMOVER = list(a1 = A), {
                              Pr_POP_SYLAB[s1,y1,lc1,a1,b2] * Pr_A_SYL[a1,s1,y1,lc1]
                            })  
    
    Pr_POP_SYLSYL <- autoloop(
      s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, s2=SEXES, y2=SAMPY, lc2=LENGTH_CLASSES,
      SUMOVER=list(a1= A, a2= A), 
    {
      b2 <- y2 - a2
      (b2 >= POPY[1]) *
      Pr_POP_SYLAB[ s1, y1, lc1, a1, b2 |> clamp(POPY)] *
      Pr_A_SYL[ a1, s1, y1, lc1] *
      Pr_A_SYL[ a2, s2, y2, lc2]
    })
    
    E_POP_SYLSYL <- n_comps_POP_SYLSYL * Pr_POP_SYLSYL 
    
    nll <- -sum(dpois(c(n_POP_SYLSYL),c(E_POP_SYLSYL),log = T),na.rm = T)
    
  ##Array for average fecundity by sex and age
  #fecun = array(0,c(50,2))
  #for()
    #!#!#!# add this to the big model once we make the age data
    ## Post-hoc Age addition for POPs and TKPs
    for (j in 1:nrow(POPage)) {
      # call SYLSYL group and age for both pairs
      age_s1 <- agePairs$s1[i]
      age_y1 <- agePairs$y1[i]
      age_l1 <- agePairs$l1[i]
      age_s2 <- agePairs$s2[i]
      age_y2 <- agePairs$y2[i]
      age_l2 <- agePairs$l2[i]
      
      age_a1 <- agePairs$a1[i]
      age_a2 <- agePairs$a2[i]
      
      ## grabbing SYL probs for both ages
      Pr_A_SYL_age1 <- Pr_A_SYL[age_a1,age_s1,age_y1,age_l1]
      Pr_A_SYL_age2 <- Pr_A_SYL[age_a2,age_s2,age_y2,age_l2]
      
      # pr for both ages
      Pr_SYL_a1a2 <- Pr_A_SYL_age1 * Pr_A_SYL_age2
      ## grabbing relevant SYLSYL prob
      Pr_POP_SYLSYL_tmp <- Pr_POP_SYLSYL1[age_s1,age_y1,age_l1,age_s2,age_y2,age_l2]
      
      age_b2 <- age_y2 - age_a2
      nll <- nll - (Pr_POP_SYLAB[age_s1, age_y1, age_l1, age_a1, age_b2] * Pr_SYL_a1a2)/Pr_POP_SYLSYL_tmp
    }
  
  REPORTO( N, Z, E_POP_SYLSYL, fec_sa)
  ##Return nll
  nll
})

##FOR NOW DISABLE recruitment parameters
new_f(parm)

tmbmap = list(log_rec=as.factor(rep(NA,length(parm$log_rec))))

testo = MakeADFun(new_f,parm,random=c("log_rec"),map = tmbmap)

testo$gr()
repp <- testo$report()
opt <- nlminb(testo$par, testo$fn, testo$gr)

sdr <- sdreport(testo)

pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")

## create data object with observed data
dat <- list()
# simulation census size
dat$Nobs <- CS_180_195
dat$A <- nrow(unique(dat$Nobs[,3]))
dat$Y <- nrow(unique(dat$Nobs[,1]))

parm <- list()
parm$LogN <- numeric(nrow(dat$Nobs[,1]))
parm$LogN_yr0 <- numeric(nrow(unique(dat$Nobs[,1])))
parm$logsd <- numeric(nrow(unique(dat$Nobs[,3])))
parm$meanS_F 

f <- function(parm) {
  getAll(dat,parm)
  
  #create other necessary data and blank objects
  count <- nrow(Nobs)
  sd <- exp(logsd)
  nll <- 0
  LogN <- log(LogN)
  LogN_yr0 <- log(LogN_yr0)
  
  for (i in 1:count) {
    #index year and age for the observation
    y <- as.numeric(Nobs[i,1] - min(Nobs[,1]) + 1)
    a <- as.numeric(Nobs[i,3])
    
    #for recruitment, set as rw?? This isn't right
    # have a separate 'recruit' parameter (LogN_yr0), should this be one value?
    if(a == 1){
      pred <- LogN_yr0[y]
      LogN[i] <- pred
    }
    #for age > 1, adjust based on average survival for m/f
    else{
      #!# need to calculate average survival first and then adjust
      pred <- LogN[i-1]
    }
    
    # calculate average survival based on LogN
    # also calculate fecundity here
    
    # adjust observations based on survival and fecundity
    for (i in 1:count) {
      
      
      # stand-in likelihood comparing observation to real
      nll <- nll - dnorm(log(as.numeric(Nobs[i,5])), exp(pred), sd[a], log=TRUE)
    }
    
    
  }

  REPORT(LogN)
  
  return(nll)
}

obj <- MakeADFun(f,parm, random="LogN", silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")






