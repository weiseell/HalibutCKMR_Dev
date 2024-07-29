### Just simulated Data input pop dyn test
# random-walk abundance
library(RTMB)

#load simulated data inputs
load("RTMB_Input/CensusSize_Sim_071224.rda")
load("RTMB_Input/POPs_df_072224.rda")

#source functions for RTMB model
source("ModelFunction/prob_la.R")

Narray = array(CS_181_195$abundance,c(50,15,2))
dat <- list()
dat$Narray <- Narray
dat$A <- 30
dat$Y <- 15
dat$Lvec <- seq(10,220,10)
dat$la_key <- read.table("RTMB_Input/LA_MeanSD.txt",header = T, sep = "\t",stringsAsFactors = F)

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
parm$log_fecF <- c(0.1,0.1)
parm$log_fecM <- c(0.1,0.1)

parm$log_fecF <- c(log(0.019),log(3.624))
parm$log_fecM <- c(log(0.019),log(3.624))

make_fecundity <- function(alpha,bexp,len){
  alpha*len^bexp
}

new_f <- function(parm){
  getAll(dat,parm)

  #exp the sd
  rec_sd = exp(log_rec_sd)
  ##An array for male/female fecundity
  log_fec = parm$log_Z
  fecM <- exp(log_fecM)
  fecF <- exp(log_fecF)
  ##Array for average fecundity by sex and age
  ## create prob_by_length array and then fecun array
  fecun = array(0,c(A,2))
  prob_len_at_age <- prob_la(Lmax = max(Lvec),
                             A = A, binsize = 10, 
                             dat = la_key)
  
  for(s in 1:2){
    if(s == 1){
      alpha = fecM[1]
      bexp = fecM[2]
    }
    if(s == 2){
      alpha = fecF[1]
      bexp = fecF[2]
    }
    fec_denom <- 0
    for (a in 2:A) {
      runfec <- 0
      for (i in 1:length(Lvec)) {
        #set L from length vector
        l <- Lvec[i]
        runfec <- runfec + prob_len_at_age[a-1,i,s] * make_fecundity(alpha,bexp,l)
        fec_denom <- fec_denom + runfec
      }
      fecun[a,s] <- runfec
    }
    fecun[,s] <- fecun[,s]/fec_denom
  }
  
  
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
  
  ## nll using CKMR data
  Pr_MOP_SYLASYA <- array(0,c(2,Y,length(Lvec),A,2,Y,A))
  # generating the POP probabilities
  ## issue: creation of negative birth years
  for (s1 in 1:2) { # Sex of Parent
    for (y1 in 1:Y) { # Sampling Years for Parents
      for (l1 in 1:length(Lvec)) { # Length of Parent
        for (a1 in 1:A) { # Age of Parent
          for (s2 in 1:2) { #sex of the offspring
            for (y2 in 1:Y) { # Sampling Year for Offspring
              for (a2 in 1:A) { # Age of Offspring
                #age of offspring
                B2 <- y2-a2
                B1 <- y1-a1
                #probability of zero scenarios
                if(y1 < B2){
                  Prob <- 0
                }
                else if(B1 >= B2 | B2 <= 0) {
                  Prob <- 0
                }
                
                else{
                  #####LENGTH-BASED NOT IMPLEMENTED YET - HOW TO??##
                  #how to recalculate Linf based on a single value???
                  #L_inf_A <- l1/(1-exp())
                  #calculate the age at juvenile birth year
                  #a_b2 <- a1 = (y1-B2)
                  #calculate length at juvenile birth year
                  #how to do this???
                  #l_b2
                  #get probability for the pair
                  #Prob <- fecun(l_b2,s1)/TROgen[s1,B2]
                  
                  #Prob <- fecun[B2,s1]/TROgen[B2,s1]
                  #If parent was sampled before juv BY, need to adjust for survivas1l
                  if(y1 < B2){
                    ##Leading fraction is parents that survived to birth year
                    ##Double check my indices!
                    Prob <- (N[a1,B2,s1]/N[B1+y1,y1,s1])*fecun[a1,s1]/N[1,B2,s2]
                  }else{
                  Prob <- fecun[a1,s1]/N[1,B2,s2]
                  }
                }
              Pr_MOP_SYLASYA[s1,y1,l1,a1,s2,y2,a2] <- Prob
              }
            }
          }
        }
        #calculating the prob of being a particular age
        #!# this is for half-sibs so it is turned off for now
        #given sex, year, and length
        #Pr_A_SYL <- array(0,c(A,length(Lvec),2,Y))
        #Prob_denom <- 0
        #for (a1 in 2:A) {
          #!# again how to make this 'prob given length object??'
        #  Prob_num <- prob_len_at_age[a1-1,l1,s1]*N[a1,y1,s1]
        #  Prob_denom <- Prob_denom + Prob_num
        #  Pr_A_SYL[a1,l1,s1,y1] <- Prob_num
        #}
        
        #for (a1 in 2:A) {
          # what does this do????
        #  Pr_A_SYL[a1-1,l1,s1,y1] = Pr_A_SYL[a1-1,l1,s1,y1]/Prob_denom
        #}
      }
    }
  }
  
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






