### Just simulated Data input pop dyn test
# random-walk abundance
library(RTMB)
library(offarray)

#load simulated data inputs
load("Inputs/CensusSize_Sim_071224.rda")
load("Inputs/POPs_df_072224.rda")

#source functions for RTMB model
source("ModelFunction/prob_la.R")

Narray = array(CS_181_195$abundance,c(50,15,2))
dat <- list()
dat$Narray <- Narray
dat$A <- 2:30
dat$POPY <- 1:15
dat$SAMPY <- 9:15
dat$Lvec <- seq(10,220,10)
dat$LENGTH_CLASSES <- 1:length(dat$Lvec)
dat$SEXES <- 1:2
dat$la_key <- read.table("Inputs/LA_MeanSD.txt",header = T, sep = "\t",stringsAsFactors = F)

#number of POPs and comps as array (from script 2 that processes sim results)
dat$n_POP_sylasya <- nonzero

## create parameter list
parm <- list()
#log recruitement - random effect length Y
parm$log_rec <- log(Narray[1,,1] + Narray[1,,2])
parm$log_rec_sd = log(0.3)
#parm$log_rec <- rep_len(x = 0.1, length.out = dat$Y)
#total log abundance - split into Y/A abundance within the model
#parm$log_init_abundance <- 0.1
parm$log_init_abundance <- log(sum(Narray[-1,1,]))
# average log survival for M and F
parm$log_Z <- c(0.1,0.1)
#a and b values for the fecundity equations
#!# fix these for now
parm$log_bexp <- c(0.1,0.1)

make_fecundity <- function(s,len){
  (len/150)^bexp[s]
}

new_f <- function(parm){
  getAll(dat,parm)

  #exp the sd
  rec_sd = exp(log_rec_sd)
  ##An array for male/female fecundity
  bexp <- exp(log_bexp)
  
  ##Array for average fecundity by sex and age
  ## create prob_by_length array and then fecun array
  fecun = array(0,c(2,A))
  prob_len_at_age <- prob_la(Lmax = max(Lvec),
                             A = A, binsize = 10, 
                             dat = la_key)
  
  fec_sa <- autoloop(
      s=SEXES, a=A, 
      SUMOVER=list( lc=LENGTH_CLASSES), {
    l <- Lvec[lc]
    make_fecundity(s, l) * prob_len_at_age[l,a,s]
  })
  
  ##Plus group for later
  N = array(0,c(A,Y,2))

  nll = 0
  ##Random walk for recruitment
  #for(y in 2:Y){
  #  nll = nll - dnorm(log_rec[y],log_rec[y-1],rec_sd,TRUE)
  #}
  
  ##Put in recruitment
  for(s in 1:2){
    N[1,,s] = exp(log_rec)/2
  }
  
  Z = array(0,c(A,Y,2))
  ##Male Mortality
  Z[,,1] = exp(log_Z[1])
  ##Female
  Z[,,2] = exp(log_Z[2])
  
  cumsurv_m = rev(cumsum(Z[-1,1,1]))
  cumsurv_f = rev(cumsum(Z[-1,1,2]))
  total_surv = sum(cumsurv_m+cumsurv_f)
  
  N[-1,1,1] = (cumsurv_m/total_surv)*exp(log_init_abundance)
  N[-1,1,2] = (cumsurv_f/total_surv)*exp(log_init_abundance)
  
  
  for (y in 2:Y) {
    for (a in 2:A) {
      for (s in 1:2) {
        N[a,y,s] <- N[a-1,y-1,s]*exp(-Z[a-1,y-1,s])
      }
    }
  }
  
  TRO_SY <- autoloop(
      s=SEXES, y=POPDYN_YEARS,
      SUMOVER( a=AGES),{
    N[a,y,s] * fec_sa(s,a)    
  })
  
  inv_TRO_SY <- 1/TRO_SY
  
  ## nll using CKMR data
  Pr_MOP_SYLASYA <- array(0,c(2,Y,length(Lvec),A,2,Y,A))
  # generating the POP probabilities
                # This does all ages for IDEAl measurement
                # autoloop() both CREATES the array and FILLS IT IN
                # looping over all the "equals" at the start
                # This does the plus-group wrong. We will fix immediately afterwards...
                Pr_POP_SYLAB <- autoloop( 
                    s1=SEXES, y1=SAMP_YEARS, lc1=LENGTH_CLASSES, a1=AGES,
                    b2=POPDYN_YEARS, {
                  #sd that parent length is currently from the mean
                  
                  l1 <- actual_l[ lc1]
                  sd1 <- (l1 - la_means_SA[s1,a1])/la_sd_SA[s1,a1]
                  
                  #age of parent when off is born
                  a1_at_B2 <- a1 - (y1 - b2)
                  #length of parent when offspring is born
                  l1_at_B2 <- la_means_SA[s1,a1_at_B2] + sd1 * la_sd_SA[s1,a1_at_B2]
                  
                  #!# switch fecundity input from age-based to length-based in the fecundity
                  Prob <- 
                    (y1 >= b2) * # otherwise Molly was dead before Dolly born
                    (a1_at_B2 >= 2) *
                    feclen_fun(s1,l1_at_B2) * inv_TRO_SB[s1,B2]
                })
                
                # PLUS GROUP FIXUP HERE
                
    num_Pr_A_SYL <- autoloop(
      a=AGES, s=SEXES, y=SAMPY, lc=LENGTH_CLASSES,
      prob_len_at_age[s,l,a] * N[s,y,a]
    )              

    denom_SYL <- sumover( 'a', num_Pr_A_SYL)
    
    Pr_a_SYL <- autoloop(
      a=AGES, s=SEXES, y=SAMPY, lc=LENGTH_CLASSES,
      num_Pr_a_SYL[ a, s, y, l] / denom_SYL[ s, y, l]
    )
    
    Pr_POP_SYLSYL <- autoloop(
      s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, s2=SEXES, y2=SAMPY, lc2=LENGTH_CLASSES,
      SUMOVER=list( a1= AGES, a2= AGES), 
    {
      b2 <- y2 - a2
      Pr_POP_SYLAB[ s1, y1, lc1, a1, b2] *
      Pr_a_SYL[ a1, s1, y1, lc1] *
      Pr_a_SYL[ a2, s2, y2, lc2]
    })
 
  ## adding likelihood based of POP composition
  # need to encompass all parts of the Prob vector
  # this one only has parent length, may need to add uncertain age?
  # also need to adjust these dimensions for what the non-zero comparisons are
  for (s1 in 1:2) {
    for (y1 in 1:Y) {
      for (l1 in 1:length(Lvec)) {
        for (a1 in 1:A) {
          for (s2 in 1:2) {
            for (y2 in 1:Y) {
              for (a2 in 1:A) {
                #these arrays are a PROBLEM, will need to check for future sims and emp data.
                #data once this version runs
                
                #running df version
                #default is n_comps/n_POPs is 0
                #check df for nonzero values
                n_comps <- 0
                n_POPs <- 0
                co <- 1
                while (co < nrow(n_POP_sylasya)) {
                  if (n_POP_sylasya$sex.x[co] == s1 &
                      n_POP_sylasya$SampYear.x[co] == y1 &
                      n_POP_sylasya$lbin[co] == l1 &
                      n_POP_sylasya$AgeAtSamp.x[co] == a1 &
                      n_POP_sylasya$sex.y[co] == s2 &
                      n_POP_sylasya$SampYear.y[co] == y2 &
                      n_POP_sylasya$AgeAtSamp.y[co] == a2) {
                    n_comps <- n_POP_sylasya$ncomp[co]
                    n_POPs <- n_POP_sylasya$POP[co]
                  }
                  co <- co + 1
                }
                
                ##code for array version
                #n_comps <- n_comps_POP_sylasya[s1,y1,l1,a1,s2,y2,a2]
                #n_POPs <- n_POP_sylasya[s1,y1,l1,a1,s2,y2,a2]
                
                ## add these to the likelihood
                # number of pops and lambda is comps weighted by the probability
                # of the pairs as calculated in the above POP prob loop
                nll <- nll - dbinom(x = n_POPs,
                                   size = n_comps,
                                   prob = Pr_MOP_SYLASYA[s1,y1,l1,a1,s2,y2,a2],T)
              }
            }
          }  
        }
      }
    }
  }
  
  ##Array for average fecundity by sex and age
  #fecun = array(0,c(50,2))
  #for()
  
  REPORT(N)
  REPORT(Z)
  ##Return nll
  nll
}

##FOR NOW DISABLE recruitment parameters
tmbmap = list(log_rec=as.factor(rep(NA,length(parm$log_rec))))

testo = MakeADFun(new_f,parm,random=c("log_rec"))
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






