### Just simulated Data input pop dyn test
# random-walk abundance
library(tidyverse)
library(RTMB)
library(offarray)
library(mvbutils)
library(offartmb)

#load simulated data inputs
load("Inputs/CensusSize_Sim_090325.rda")

# read in otolith data
lengthage <- read.csv("E:/ClonedRepositories/ForTheHalibut/Input/Armsworthy&Campana_ladata_nofilter.csv")

#source functions for RTMB model
source("ModelFunction/prob_la.R")
source("ModelFunction/growthFunGen.R")
source("ModelFunction/LengthBinInputs.R")

##read in the sim check
simCheck =  readRDS("simCheck.rds")

# read in the simulated POP and HSP inputs
rel_counts = read.csv("./Inputs/Relate_ncomp_df_090425.csv")

####
## dat inputs for the model
####
Narray = array(CS_181_195$abundance,c(50,15,2))
dat <- list()
dat$Narray <- Narray
dat$A <- 2:30
dat$POPY <- 1:15
dat$SAMPY <- 10:15
dat$SEXES <- 1:2
dat$Lvec <- 1:10
dat$male <- 1
dat$female <- 2
dat$PopYear1 <- dat$POPY[1]
dat$la_key <- lengthage %>% 
  filter(OtoAge <= max(A)) %>% 
  group_by(ObsSex,OtoAge) %>% 
  summarise(Means=mean(Length),SDs=sd(Length)) %>% 
  rename(sex=ObsSex,AgeClass=OtoAge,mean=Means,sd=SDs)

###CREATE POP/HSP offarrays
n_POP_SYLAB <- offarray(x=0,dimseq= list(s1=dat$SEXES,y1=dat$SAMPY,lc1=dat$Lvec,a1=dat$A,b2=dat$POPY))

n_comps_POP_SYLAB <- offarray(x=0,dimseq= list(s1=dat$SEXES,y1=dat$SAMPY,lc1=dat$Lvec,b2=dat$POPY))

N_TKP_BB <- offarray(x=0,dimseq= list(b1=dat$POPY,b2=dat$POPY))

N_comps_TKP_BB <- offarray(x=0,dimseq= list(s1=dat$SEXES,y1=dat$SAMPY,lc1=dat$LENGTH_CLASSES,s2=dat$SEXES,y2=dat$SAMPY,lc2=dat$LENGTH_CLASSES))

# function to load data into offarrays
fill_offcrud <- function(count_df,count_array,target_col="POP"){
    count_df = count_df[count_df[,target_col] != 0,]
    for(i in 1:nrow(count_df)){
        count_array[count_df$sex.x[i],count_df$SampYear.x[i],count_df$lbin.x[i],count_df$sex.y[i],count_df$SampYear.y[i],count_df$lbin.y[i]] = count_df[i,target_col]
    }
    count_array
}

#add offarrays to data input
dat$n_POP_SYLSYL = fill_offcrud(rel_counts,n_POP_SYLSYL,"POP")
dat$n_comps_POP_SYLSYL = fill_offcrud(rel_counts,n_comps_POP_SYLSYL,"ncomp")
dat$N_TKP_SYLSYL = fill_offcrud(rel_counts,N_TKP_SYLSYL,"SOPs")


## SIMULATED HAPLOTYPE FREQS
## define haplotype structure
nHap <- 7
HapFreq <- 1:nHap
HapFreq <- HapFreq/sum(HapFreq)
dat$HapFreq <- HapFreq

#split la_key into arrays of means and sds
raw_la_means_SA <- dat$la_key %>% 
  arrange(sex,AgeClass) %>% 
  dplyr::select(sex,AgeClass,mean) %>% 
  filter(AgeClass <= 30) %>% 
  spread(AgeClass,mean)
raw_la_means_SA <- raw_la_means_SA[,-1]

raw_la_sd_SA <- dat$la_key %>% 
  arrange(sex,AgeClass) %>% 
  dplyr::select(sex,AgeClass,sd) %>% 
  filter(AgeClass <= 30) %>% 
  spread(AgeClass,sd)
raw_la_sd_SA <- raw_la_sd_SA[,-1]

dat$la_means_SA <- offarray(as.matrix(raw_la_means_SA),dimseq=list(s=dat$SEXES, a=dat$A))
dat$la_sd_SA <- offarray(as.matrix(raw_la_sd_SA),dimseq=list(s=dat$SEXES, a=dat$A))

# create prob len at age array
#!# need to figure out what to do about how to reconvert this to quantile lengths now
raw_prob_len_at_age <- prob_la(Lmax = max(dat$Lvec),
                               ages = dat$A,
                               binsize = 10,
                               dat = dat$la_key)

# Make it into offarray
dat$prob_len_at_age <- offarray(raw_prob_len_at_age,
                                dimseq=list( s=dat$SEXES, lc=dat$Lvec, a=dat$A))

### Fecundity function fix
#defining shape and scale parameters for pgamma in the
#fecundity functions
dat <- within(dat, {
  la_shape_SA <- (la_means_SA/la_sd_SA)^2
  la_scale_SA <- la_sd_SA^2/la_means_SA
})

####
## create parameter list
####
parm <- list()
#log recruitment - random effect length popdyn years
parm$rw_log_rec <- rep(log(10000),length(Narray[1,,1]))
parm$rw_log_rec_sd = log(0.3)
parm$noise_logrec_dev <- Narray[1,,1] * 0
parm$noise_rec_dev_sd_log = log(0.3)

parm$log_rec <- c(log(10000),log(10000))
#total log abundance - split into Y/A abundance within the model
#currently have the summed abundance of year 181 in the simulation
parm$log_init_abundance <- log(74884)
# average log survival for M and F
parm$log_Z <- log(c(0.18,0.18))
#a and b values for the fecundity equations
#grabbed the bexp value from the simulation
#currently there is no a value in the model
parm$log_bexp <- log(c(3.624,3.624))

#no lucky litter analog in the simulation
parm$log_lucky_litter_par <- 2

#geometric slope of recruitment
##parm$logslopea0 <- log(0.2)
#sex ratio - sim has equal sex ratio
parm$ppn_female <- 0.5

#!# Start here with loading the environment!!!

## function input to RTMB
new_f <- function(parm) reclasso(by=parm, {
  getAll(dat,parm) 
  
  len_at_age1 = splinefun(x=pdat1$AgeClass,y=pdat1$pred)
  len_at_age2 = splinefun(x=pdat2$AgeClass,y=pdat2$pred)
  
  
  sd_at_len1 = splinefun(x=sdpredr1$mean,y=sdpredr2$sd)
  sd_at_len2 = splinefun(x=sdpredr1$mean,y=sdpredr2$sd)
  
  
  #exp logged parameter values
  rw_rec_sd = exp(rw_log_rec_sd)
  bexp <- exp(log_bexp)
  ##rec <- exp(log_rec)
  noise_rec_dev_sd <- exp(noise_rec_dev_sd_log)
  lucky_litter_par <- exp(log_lucky_litter_par)
  slopea0 <- mean(exp(log_Z))
  ##slopea0 <- exp(logslopea0)
  
  ##Array for average fecundity by sex and age
  # centered around 150cm
  make_fecundity <- function(s,len){
    (len/150)^bexp[s]
  }
   
  ## create fecundity arrays
  fec_sa <- autoloop(
    s=SEXES, a=A, 
    SUMOVER=list(lc=Lvec), {
      make_fecundity(s,LengthQ_SAL[SLICE=s,SLICE=a,SLICE=lc]) * prob_len_at_age[s,lc,a]
    })
  
  fec_sa_quant <- autoloop(
    s=SEXES, a=A, qq=(1:nquant), {
      #calculating length for quantile at mean/sd
      l_q <- qgamma(qq/(nquant + 1),
                    shape = la_shape_SA[s,a],
                    scale = la_scale_SA[s,a])
      
      #make fecundity for length given prob
      make_fecundity(s,l_q) * (1/nquant)
    })

  # initialize total abundance array
  N = offarray(0,dimseq = list(SEXES=SEXES,POPY=POPY,AGE=A))
  
  # initialize nll
  nll = 0
  
  # random walk with noise for recruitment
  ##nll <- nll - sum(dnorm(diff(rw_log_rec), mean=0, sd=rw_rec_sd))
  ##nll <- nll - sum(dnorm(noise_logrec_dev, mean=0, sd=noise_rec_dev_sd))
  for(i in 2:length(rw_log_rec)){
    nll <- nll - dnorm(rw_log_rec[i],rw_log_rec[i-1],rw_rec_sd)
  }
  
  
  # cumulating recruitment
  ##cumul_rw_logrec <- cumsum(rw_log_rec)
  
  #per year change in recruitment
  ##rec_mul <- exp(cumul_rw_logrec + noise_logrec_dev) 
  ##rec_mul <- exp(cumul_rw_logrec)
    
  
  ##Put in recruitment into N
  for(s in 1:2){
    #N[s,,2] = rec[s] * rec_mul
    N[s,,2] = exp(rw_log_rec)
  }
  
  # initizalize mortality vector
  Z = 0 * N
  
  ##Male Mortality
  Z[1,,] = exp(log_Z[1])
  ##Female
  Z[2,,] = exp(log_Z[2])
  
  ## calculating the first year for all ages
  # initialize
  AAA = A[-1]
  Ny0 <- offarray(0, dimseq = list(AGE=AAA))
  Ny0[min(AAA)] <- 1
  
  # add the rest of the years using slope a0
  for (a in (min(AAA)+1):max(AAA)) {
    Ny0[a] <- Ny0[a-1] * exp(-slopea0)
  }
  
  # adjust the plus group
  Ny0[max(AAA)] <- Ny0[max(AAA)]/(1-slopea0)
  Ny0 = Ny0/sum(Ny0)
  
  # adjust the whole vector by total initial abundance
  Ny0 <- Ny0 * exp(log_init_abundance)
  
  # load vector into N
  # adjust male/female based on sex ratio parameter
  for(aa in (min(AAA)):max(AAA)){
  N[female,min(POPY),aa] <- Ny0[aa] * ppn_female
  N[male,min(POPY),aa] <- Ny0[aa] * (1-ppn_female)
  }
  ## initialize mean plus group age
  Abar_plus <- offarray(0,dimseq = list(SEXES=SEXES,POPY=POPY))
  Abar_plus[,1] <- max(A) + 1/(1-exp(-Z[,SLICE=1,SLICE=max(A)]))
  
  #set adults as years without recruitment
  adults <- A[-1]
  
  # loop to fill in the rest of N
  for (y in 2:length(POPY)) {
    ## run abundance for the next year based on mortality
    N[,y,adults] <- autoloop(s = SEXES, a = adults,{
      N[s,y-1,a-1]*exp(-Z[s,y-1,a-1])
    })
    
    #adjust the plus group
    N[,y,max(A)] <- N[,SLICE=y,SLICE=max(A)] + N[,SLICE=y-1,SLICE=max(A)] * exp(-Z[,SLICE=y-1,SLICE=max(A)])
    
    #add the Abar for this year based on the plus group
    Abar_plus[,y] <- c(((Abar_plus[,SLICE=y-1]+1) * N[,SLICE=y-1,SLICE=max(A)] + 
                          max(A) * N[,SLICE=y-1,SLICE=max(A)-1])/
                         (N[,SLICE=y-1,SLICE=max(A)] + N[,SLICE=y-1,SLICE=max(A)-1]))
    
  }
  
  ## calculate survival array to deal with plus group in probabilities
  ## used in half-sib probability
  # make a new max age vector for plus group
  maxage <- max(A) + max(POPY)
  superAge <- 2:maxage
  
  ## run abundance for plus group years
  # make extra years offarray
  superN <- offarray(0,dimseq = list(SEXES=SEXES,POPY=POPY,AGE=superAge))
  superN[,,A] <- N
  superN[,2:max(POPY),max(A)] <- 0
  # add N for those extra years
  for (y in 2:length(POPY)) {
    # needed to vectorize Z and superN slices to multiply them
    superN[,SLICE=y,max(A):max(maxage)] <- c(superN[,SLICE=(y-1),(max(A):max(maxage))-1]) * 
      c(exp(-Z[,SLICE=(y-1),SLICE=max(A)]))
  }
  
  ## calculate survival past 30 with superN
  Pr_Surv_SYAY <- autoloop(
    spp=SEXES,y1=POPY,app=superAge,y2=POPY,{
      superN[spp,y2,(app + (y2-y1)) |> clamp(superAge)]/superN[spp,y1,app]
    }
  )
  
  #calculate reproductive output per sex/age group
  TRO_SY <- autoloop(
    s=SEXES, y=POPY,
    SUMOVER=list(a=A),{
      N[s,y,a] * fec_sa[s,a]    
    })
  
  #making the reciprocal because RTMB hates dividing stuff
  recip_TRO_SY <- 1/TRO_SY
  ####
  ## POP nll using CKMR data ###
  ####
  # generating the POP probabilities
  # This does all ages for IDEAL measurement
  # autoloop() both CREATES the array and FILLS IT IN
  # looping over all the "equals" at the start
  # This does the plus-group wrong. We will fix in the next section
  Pr_POP_SYLAB <- autoloop( 
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, a1=A,
    b2=POPY, {
      # generate probability given fecundity of par
      # and TRO at the offspring birth year
      Prob <- (y1 >= b2) * # otherwise Molly was dead before Dolly born
        (a1_at_B2 >= 2) *
        make_fecundity(s1,LengthQ_SYLAB[s1,y1,lc1,a1,b2]) * recip_TRO_SY[s1,b2]
      
      Prob
    })
  
  # PLUS GROUP FIXUP HERE
  # adjusting probabilities to account for plus group
  # age of parent is 30+ instead of only 30
  #Pr_POP_SYLAB_plus <- autoloop(s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES,
  #                              b2=POPY, SUMOVER=list(q=1:nquant), {
                                  ##sd that parent length is currently from the mean
                                  # calculate what a1 is based on the length
  #                                a1 <- qexp(q/(nquant + 1)) * (Abar_plus[s1,y1] - max(A)) + max(A)
                                  
  #                                l1 <- Lvec[lc1]
                                  #!# CAN FIX THIS IN THE FUTURE BY CHANGING THE MEANS/SD FOR LENGTH/AGE
                                  # TO A FUNCTION BASED ON LA DATA
  #                                sd1 <- (l1 - len_at_age1(a1))/sd_at_len1(a1)
                                  
                                  #age of parent when off is born
  #                                a1_at_B2 <- a1 - (y1 - b2)
  #                                l1_at_B2 <- len_at_age1(a1_at_B2) + 
  #                                  sd1 * sd_at_len1(a1_at_B2)
                                  
                                  #probability based on par fecundity
                                  # and TRO at off birth
  #                                Prob <- (y1 >= b2) * # otherwise Molly was dead before Dolly born
  #                                  (a1_at_B2 >= 2) * 
  #                                  (l1_at_B2 > 0) * # otherwise length shows that Molly wouldn't have been born yet
  #                                  make_fecundity(s1,l1_at_B2) * recip_TRO_SY[s1,b2]
                                  
  #                                Prob
  #                              })
  
  # switching plus group prob
  #Pr_POP_SYLAB[,,,max(A),] <- Pr_POP_SYLAB_plus
  
  
  #expected number of POPs given probability and number of comparisons
  E_POP_SYLAB <- n_comps_POP_SYLAB * Pr_POP_SYLAB
  
  #summing all the calculated likelihoods and adding them to the total
  nll <- nll - sum(dpois(c(n_POP_SYLAB),c(E_POP_SYLAB),log = T),na.rm = T)
  
  ##REPORTO(N, Z, E_POP_SYLSYL, E_N_TKP_SYLSYL, fec_sa)
  REPORT(N)
  REPORT(Z)
  REPORT(Pr_POP_SYLAB)
  REPORT(TRO_SY)
  REPORT(Pr_Surv_SYAY)
  REPORT(Pr_TKP_SYLSYL)
  
  ##Return nll
  nll
})

##FOR NOW DISABLE recruitment parameters
tmpout <- new_f(parm)

tmbmap = list(noise_logrec_dev=as.factor(rep(NA,length(parm$noise_logrec_dev))),
              noise_rec_dev_sd_log = as.factor(NA),
              ppn_female = as.factor(NA),
              log_lucky_litter_par = as.factor(1),
              log_Z = as.factor(c(1,1)),
              log_bexp = as.factor(c(NA,NA)))

testo = MakeADFun(new_f,parm,random=c("rw_log_rec"),silent=FALSE,map=tmbmap)

#badpar <- c(-1.20352,10.0034,9.99660,11.0040,-1.59915,-1.58924,0.00239148,0.0397631,0.0269279,-1.63183,0.418786)
#testo$fn(badpar)
#testo$gr(badpar)

##grady = testo$gr()
##names(grady) = names(testo$par)


repp <- testo$report()

Pr_POP_SYLAB = repp$Pr_POP_SYLAB
Pr_POP_SYAB = apply(Pr_POP_SYLAB,c(1,2,4,5),sum)



dadC = simCheck$dadCdf |>
    mutate(pdy = as.numeric(as.character(pdy))) |>
    filter(pdy >= 190, pdy <= 195) |>
    filter(kby >= 190, kby <= 195) |>
    mutate(pdy2 = factor(pdy,levels=190:195,labels=10:15),kby2=factor(kby,levels=190:195,labels=10:15))

momC = simCheck$momCdf |>
    mutate(mdy = as.numeric(as.character(mdy))) |>
    filter(mdy >= 190, mdy <= 195) |>
    filter(kby >= 190, kby <= 195) |>
    mutate(mdy2 = factor(mdy,levels=190:195,labels=10:15),kby2=factor(kby,levels=190:195,labels=10:15))

Pr_POP_SYABdf = as.data.frame.table(Pr_POP_SYAB)

Pr_POP_YBdf = split(Pr_POP_SYBdf,Pr_POP_SYBdf$s1)



opt <- nlminb(testo$par, testo$fn, testo$gr,control = list(trace=0,iter.max=1000,eval.max=1000))

expsdr <- sdreport(testo)

pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")

