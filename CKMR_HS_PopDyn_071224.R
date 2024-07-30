### Just simulated Data input pop dyn test
# random-walk abundance

#load libraries
library(RTMB)

setwd("~/Desktop/HalibutCKMR_Dev")

#load simulated data inputs
load("Inputs/CensusSize_Sim_071224.rda")
load("Inputs/Ncomp_HSs_yaya_071724.rda")
load("Inputs/NHSs_syaya_071724.rda")

#source functions for RTMB model
source("ModelFunction/prob_la.R")
source("ModelFunction/quant_la.R")

Narray = array(CS_181_195$abundance,c(50,15,2))
dat <- list()
dat$Narray <- Narray
dat$A <- 30
dat$Y <- 15
dat$Lvec <- seq(10,220,10)
dat$nquant <- 4
dat$la_key <- read.table("Inputs/LA_MeanSD.txt",header = T, sep = "\t",stringsAsFactors = F)

#number of POPs and comps as array (from script 2 that processes sim results)
dat$n_comps_HS_yaya <- n_comps_HS_yaya
dat$n_SO_kins_syaya<- n_SO_kins_syaya

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
  for(y in 2:Y){
    nll = nll - dnorm(log_rec[y],log_rec[y-1],rec_sd,TRUE)
  }
  
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
  ## Addition of information from HSPs
  # I have NO IDEA how this is working
  # basically I just transcribed this bit from the outline
  Pr_OHSP_sbb <- array(0,c(2,Y,Y))
  quant_len_at_age <- quant_la(A,nquant,la_key)
  for (sp in 1:2) {
    if(sp == 1){
      alpha = fecM[1]
      bexp = fecM[2]
    }
    if(sp == 2){
      alpha = fecF[1]
      bexp = fecF[2]
    }
    
    for (B1 in 1:Y) {
      for (B2 in 1:Y) {
        ##!## need to make sure all probs are less than 1
        ## have a fecun issue in the quantiles section I thint
        ## not currently sure how to fix
        Pr_HSP_running <- 0
        ## eliminate same-cohort half-sibs? How?
        if (B1 == B2) {
          Pr_HSP_running <- 0
        } ## if B1 is younger than B2, prob is zero
        ## don't add to running prob
        if(B1 > B2){
          Pr_HSP_running <- 0
          
        }else{
          for (app in 2:A) {
            ## if parent had to survive past max age
            ## don't add to running prob
             if(app + (B2-B1) > A){
              Pr_HSP_running <- Pr_HSP_running
              } else {
              ## prob of parent surviving from B1 to B2
              Pr_surv <- N[app + (B2-B1), B2, sp]/N[app,B1,sp]
              
              pp_ERO_B2 <- 0
              denom <- 0
              
              for (qu in 1:nquant) {
                #!# what is num_qs and the quantile thing?????
                # what is length_at_age_q?? I think it is an array but not sure how to make it
                #right now it's written as a data frame that is age x length, would need to be sex-separated
                #or could be a 3D array 
                
                Pr_pp_is_qth_q_1 <- 1/(nquant*make_fecundity(alpha,bexp,quant_len_at_age[app-1,qu,sp]))
                denom <- denom + Pr_pp_is_qth_q_1
                
                ##ERO for B2
                
                pp_ERO_B2_if_q <- make_fecundity(alpha,bexp,quant_len_at_age[app+B2-B1-1,qu,sp])
                #ERO for all potential birth years
                #weighted by prob of length at age?
                pp_ERO_B2 <- pp_ERO_B2 + Pr_pp_is_qth_q_1 * pp_ERO_B2_if_q
              }
              # dividing ERO by Pr - denominator across all potential quantiles
              pp_ERO_B2 <- pp_ERO_B2/denom
              
              # probability of shared parent??? Not exactly sure how all this comes together
              #!# not right N?
              Pr_pp_is_1s_mum <- fecun[app,sp]/N[1,B1,sp]
              
              #putting the HSP probability together
              Pr_HSP_running <- Pr_HSP_running + N[app,B1,sp] * Pr_pp_is_1s_mum * Pr_surv * (pp_ERO_B2/N[1,B1,sp])
              
              ## compensate for same-cohort HSs (lucky litter)
              #need to make a vector of some kind called lucky_litter_factor
              #Pr_HSP_running <- Pr_HSP_running * lucky_litter_factor(sp)
            }
          }
        }
        Pr_OHSP_sbb[sp,B1,B2] <- Pr_HSP_running
      }
    }
  }
  
  ## Add GSP cases
  ## Add mtDNA Data as a post-hoc probability
  
  ## Likelihood with age composition data
  #!# need an array input with the number of otoliths at each given age/year/length/sex
  #!# also need array that is not separated by age, just length
  #for (s in 1:2) {
  #  for (y in 1:Y) {
  #    for (l in 1:L_max) {
  #      for (a in 1:A) {
          # neg log-likelihood adjustment for oto ages
  #        nll = nll - dpois(n_otos_at_age[a,l,s,y],
  #                          #total number of otoliths, adjust CKMR Prob based on that 
  #                          n_otos[l,s,y] * Pr_A_SYL[s,y,l],T)
  #      }
  #    }
  #  }
  #}
  
  ## adding likelihood based on HS composition
  for (y1 in 1:Y) {
    for (a1 in 1:A) {
      for (y2 in 1:Y) {
        for (a2 in 1:A) {
          for (sp in 1:2) {
            #!# need this array as well, number of comparisons for HSs for each y/a class for 
            #both pairs
            n_comps <- n_comps_HS_yaya[y1,a1,y2,a2]
            # number of second-order KP - half-sibs?
            #!# this array we'd know both ages in the pair & pat/mat from mtDNA
            n_SO_kins <- n_SO_kins_syaya[sp,y1,a1,y2,a2]
            #!# also need to grab the prob from the above calculation
            Pr_SO_kins <- Pr_OHSP_sbb[sp,y1,y2]
            
            # add this set to the likelihood
            nll <- nll - dbinom(n_SO_kins,n_comps,Pr_SO_kins,T)
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
tmbmap = list(log_rec=as.factor(rep(NA,length(parm$log_rec))),
              log_fecM = as.factor(c(NA,NA)),
              log_fecF = as.factor(c(NA,NA)))

testo = MakeADFun(new_f,parm,random=c("log_rec"),map = tmbmap)

repp <- testo$report()
