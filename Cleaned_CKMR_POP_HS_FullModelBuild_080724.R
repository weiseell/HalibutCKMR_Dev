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
load("Inputs/NHSs_SYLSYL_071724.rda")
load("Inputs/haps_HSs_Sim_080124.rda")
#source functions for RTMB model
source("ModelFunction/prob_la.R")
source("ModelFunction/growthFunGen.R")

####
## dat inputs for the model
####
Narray = array(CS_181_195$abundance,c(50,15,2))
dat <- list()
dat$Narray <- Narray
dat$A <- 2:30
dat$POPY <- 1:15
dat$SAMPY <- 10:15
dat$Lvec <- seq(10,220,10)
dat$LENGTH_CLASSES <- 1:length(dat$Lvec)
dat$SEXES <- 1:2
dat$nquant <- 5
dat$male <- 1
dat$female <- 2
dat$PopYear1 <- dat$POPY[1]
dat$la_key <- read.table("Inputs/LA_MeanSD.txt",header = T, sep = "\t",stringsAsFactors = F)

#number of POPs and comps as array (from script 2 that processes sim results)
dat$n_POP_SYLSYL <- n_POP_SYLSYL
dat$n_comps_POP_SYLSYL <- n_comps_POP_SYLSYL
dat$N_TKP_SYLSYL <- N_TKP_SYLSYL
dat$TKPairs <- TKPairs

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
  spread(AgeClass,mean) %>% 
  dplyr::select(-sex)

raw_la_sd_SA <- dat$la_key %>% 
  arrange(sex,AgeClass) %>% 
  dplyr::select(sex,AgeClass,sd) %>% 
  filter(AgeClass <= 30) %>% 
  spread(AgeClass,sd) %>% 
  dplyr::select(-sex)

dat$la_means_SA <- offarray(as.matrix(raw_la_means_SA),dimseq=list(s=dat$SEXES, a=dat$A))
dat$la_sd_SA <- offarray(as.matrix(raw_la_sd_SA),dimseq=list(s=dat$SEXES, a=dat$A))

# create prob len at age array
raw_prob_len_at_age <- prob_la(Lmax = max(dat$Lvec),
                               ages = dat$A,
                               binsize = 10,
                               dat = dat$la_key)

# Make it into offarray
dat$prob_len_at_age <- offarray(raw_prob_len_at_age, 
                                dimseq=list( s=dat$SEXES, lc=dat$LENGTH_CLASSES, a=dat$A))

### Fecundity function fix
#defining shape and scale parameters for pgamma in the
#fecundity functions
dat <- within( dat, {
  la_shape_SA <- (la_means_SA/la_sd_SA)^2
  la_scale_SA <- la_means_SA/la_shape_SA
})

####
## create parameter list
####
parm <- list()
#log recruitment - random effect length popdynyears
parm$rw_log_rec <- Narray[1,,1] * 0
parm$rw_log_rec_sd = log(0.3)
parm$noise_logrec_dev <- Narray[1,,1] * 0
parm$noise_rec_dev_sd_log = log(0.3)

parm$log_rec <- c(10,10)
#total log abundance - split into Y/A abundance within the model
parm$log_init_abundance <- log(sum(Narray[-1,1,]))
# average log survival for M and F
parm$log_Z <- log(c(0.2,0.2))
#a and b values for the fecundity equations
parm$log_bexp <- c(0,0)
parm$log_lucky_litter_par <- 0

#geometric slope of recruitment
parm$logslopea0 <- log(0.2)
#sex ratio
parm$ppn_female <- 0.5

## function input to RTMB
new_f <- function(parm) reclasso( by=parm, {
  getAll(dat,parm)
  
  #exp logged parameter values
  rw_rec_sd = exp(rw_log_rec_sd)
  bexp <- exp(log_bexp)
  rec <- exp(log_rec)
  noise_rec_dev_sd <- exp(noise_rec_dev_sd_log)
  lucky_litter_par <- exp(log_lucky_litter_par)
  slopea0 <- exp(logslopea0)
  
  ##Array for average fecundity by sex and age
  # centered around 150cm
  make_fecundity <- function(s,len){
    (len/150)^bexp[s]
  }
  
  ## create fecundity arrays
  fec_sa <- autoloop(
    s=SEXES, a=A, 
    SUMOVER=list(lc=LENGTH_CLASSES), {
      make_fecundity(s,Lvec[lc]) * prob_len_at_age[s,lc,a]
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
  nll <- nll - sum(dnorm(diff(rw_log_rec), mean=0, sd=rw_rec_sd))
  nll <- nll - sum(dnorm(noise_logrec_dev, mean=0, sd=noise_rec_dev_sd))
  
  # cumulating recruitment
  cumul_rw_logrec <- cumsum(rw_log_rec)
  
  #per year change in recruitment
  rec_mul <- exp(cumul_rw_logrec + noise_logrec_dev) 
  
  ##Put in recruitment into N
  for(s in 1:2){
    N[s,,2] = rec[s] * rec_mul
  }
  
  # initizalize mortality vector
  Z = 0 * N
  
  ##Male Mortality
  Z[1,,] = exp(log_Z[1])
  ##Female
  Z[2,,] = exp(log_Z[2])
  
  ## calculating the first year for all ages
  # initialize
  Ny0 <- offarray(0, dimseq = list(AGE=A))
  Ny0[min(A)] <- 1
  
  # add the rest of the years using slope a0
  for (a in (min(A)+1):max(A)) {
    Ny0[a] <- Ny0[a-1] * slopea0
  }
  
  # adjust the plus group
  Ny0[max(A)] <- Ny0[max(A)]/(1-slopea0)
  
  # adjust the whole vector by total initial abundance
  Ny0 <- Ny0 * exp(log_init_abundance)
  
  # load vector into N
  # adjust male/female based on sex ratio parameter
  N[female,min(POPY),] <- Ny0 * ppn_female
  N[male,min(POPY),] <- Ny0 * (1-ppn_female)
  
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
      # length of parent
      l1 <- Lvec[lc1]
      # quantile of parent at a1 given l1 length
      qq1 <- pgamma(l1,shape = la_shape_SA[s1,a1],scale = la_scale_SA[s1,a1])+1e-15
      qq1 <- (a1 < 5) * pgamma(l1,shape = la_shape_SA[s1,a1],scale = la_scale_SA[s1,a1])-1e-15
      qq1 <- (a1 >= 20) * pgamma(l1,shape = la_shape_SA[s1,a1],scale = la_scale_SA[s1,a1])+1e-15
      #age of parent when off is born
      a1_at_B2 <- a1 - (y1 - b2)
      
      #length of parent when offspring is born
      l1_at_B2 <- qgamma(qq1,shape = la_shape_SA[s1,a1_at_B2 |> clamp(A)],
                         scale = la_scale_SA[s1,a1_at_B2 |> clamp(A)])
      
      # generate probability given fecundity of par
      # and TRO at the offspring birth year
      Prob <- (y1 >= b2) * # otherwise Molly was dead before Dolly born
        (a1_at_B2 >= 2) *
        make_fecundity(s1,l1_at_B2) * recip_TRO_SY[s1,b2]
      
      Prob
    })
  
  # PLUS GROUP FIXUP HERE
  # adjusting probabilities to account for plus group
  # age of parent is 30+ instead of only 30
  Pr_POP_SYLAB_plus <- autoloop(s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES,
                                b2=POPY, SUMOVER=list(q=1:nquant), {
                                  ##sd that parent length is currently from the mean
                                  # calculate what a1 is based on the length
                                  a1 <- qexp(q/(nquant + 1)) * (Abar_plus[s1,y1] - max(A)) + max(A)
                                  
                                  l1 <- Lvec[lc1]
                                  #!# CAN FIX THIS IN THE FUTURE BY CHANGING THE MEANS/SD FOR LENGTH/AGE
                                  # TO A FUNCTION BASED ON LA DATA
                                  sd1 <- (l1 - len_at_age1(a1))/sd_at_len1(a1)
                                  
                                  #age of parent when off is born
                                  a1_at_B2 <- a1 - (y1 - b2)
                                  l1_at_B2 <- len_at_age1(a1_at_B2) + 
                                    sd1 * sd_at_len1(a1_at_B2)
                                  
                                  #probability based on par fecundity
                                  # and TRO at off birth
                                  Prob <- (y1 >= b2) * # otherwise Molly was dead before Dolly born
                                    (a1_at_B2 >= 2) * 
                                    (l1_at_B2 > 0) * # otherwise length shows that Molly wouldn't have been born yet
                                    make_fecundity(s1,l1_at_B2) * recip_TRO_SY[s1,b2]
                                  
                                  Prob
                                })
  
  # switching plus group prob
  Pr_POP_SYLAB[,,,max(A),] <- Pr_POP_SYLAB_plus
  
  # calculating the numerator of the Pr_A 
  num_Pr_A_SYL <- autoloop(
    a=A, s=SEXES, y=SAMPY, lc=LENGTH_CLASSES,
    prob_len_at_age[s,lc,a] * N[s,y,a]
  )              

  # generate the reciporacle denom to avoid dividing arrays
  recip_denom_SYL <- 1 / autoloop( 
    indices=dimseq( num_Pr_A_SYL)[-1],   # see num_Pr_A_SYL above for actual list
    SUMOVER=dimseq( num_Pr_A_SYL)[1], {
      num_Pr_A_SYL[a, s, y, lc]
    }   # just the "a"
   
  )
  
  # calculating the probability
  Pr_A_SYL <- autoloop(
    a=A, s=SEXES, y=SAMPY, lc=LENGTH_CLASSES, {
      num_Pr_A_SYL[a, s, y, lc] * recip_denom_SYL[s, y, lc]
    }
  )
  
  # generating pop prob summed over ages
  # using the SYL prob
  Pr_POP_SYLB <- autoloop(s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES,b2=POPY,
                          SUMOVER = list(a1 = A), {
                            Pr_POP_SYLAB[s1,y1,lc1,a1,b2] * Pr_A_SYL[a1,s1,y1,lc1]
                          })  
  
  # blurring POP prob over all age possibilities
  Pr_POP_SYLSYL1 <- autoloop(
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, s2=SEXES, y2=SAMPY, lc2=LENGTH_CLASSES,
    SUMOVER=list(a1= A, a2= A), 
    {
      b2 <- y2 - a2
      (b2 >= POPY[1]) *
        Pr_POP_SYLAB[ s1, y1, lc1, a1, b2 |> clamp(POPY)] *
        Pr_A_SYL[ a1, s1, y1, lc1] *
        Pr_A_SYL[ a2, s2, y2, lc2]
    })
  
  # flipping array for the cases when indiv 2 is the actual parent
  Pr_POP_SYLSYL2 <- autoloop(
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, s2=SEXES, y2=SAMPY, lc2=LENGTH_CLASSES, {
      Pr_POP_SYLSYL1[s2,y2,lc2,s1,y1,lc1]
    }
  )
  
  #expected number of POPs given probability and number of comparisons
  E_POP_SYLSYL <- n_comps_POP_SYLSYL * (Pr_POP_SYLSYL1 + Pr_POP_SYLSYL2)
  
  #summing all the calculated likelihoods and adding them to the total
  nll <- nll - sum(dpois(c(n_POP_SYLSYL),c(E_POP_SYLSYL),log = T),na.rm = T)
  
  ####
  ## likelihoods with Second-Order Kin Pairs ###
  ## combines grandparent/grandchild and half-sibling probabilities
  ## Half-sibling with unseen mother
  #ideal probabilities: BB
  Pr_HSP_Mat_BB <- autoloop(
    b1=POPY,b2=POPY, SUMOVER = list(app=A,qq=(1:nquant)),{
      #specify that B2 was born after B1
      (b2 >= b1) *
        ## prob that pp was father of B1
        fec_sa_quant[female,app,qq] * recip_TRO_SY[female,b1] *
        ## prob of parent surviving from B1 to B2
        Pr_Surv_SYAY[female,b1,app,b2] *
        N[male,b2,(app + (b2-b1)) |> clamp(A)]/N[male,b1,app] *
        ## prob that pp was the father of B2
        fec_sa_quant[female,(app + (b2-b1)) |> clamp(A),qq] * recip_TRO_SY[female,b2] *
        ## lucky litter factor - same cohort HS problem
        ifelse(b2==b1,
               lucky_litter_par,1)
    }
  )
  #blur into SYLSYL
  Pr_HSP_Mat_SYLSYL <- autoloop(
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, s2=SEXES, y2=SAMPY, lc2=LENGTH_CLASSES,
    SUMOVER=list(b1=POPY,b2= POPY),{
      # calculate ages based on sample and birth year
      a1 <- y1 - b1
      a2 <- y2 - b2
      
      # blur the BB probabilities into SYLSYL
      Pr_HSP_Mat_BB[b1,b2] *
        Pr_A_SYL[a1 |> clamp(A), s1, y1, lc1] * 
        Pr_A_SYL[a2 |> clamp(A), s2, y2, lc2]
    }
  )
  
  ## Half-sibling with unseen father
  #ideal probabilities: BB
  Pr_HSP_Pat_BB <- autoloop(
    balpha=POPY,bbeta=POPY, SUMOVER = list(app=A,qq=(1:nquant)),{
      ## specify this so that (hopefully) we fill in both sides of the probability triangle
      # need to do that because when we blur over length there's no guarantee that b1 < b2
      # so we need the non-zero probability in both directions
      b1 <- pmin(balpha, bbeta)
      b2 <- pmax(balpha, bbeta)
      ## prob that pp was father of B1
      fec_sa_quant[male,app,qq] * recip_TRO_SY[male,b1] *
        ## prob of parent surviving from B1 to B2
        Pr_Surv_SYAY[male,b1,app,b2] *
        N[male,b2,(app + (b2-b1)) |> clamp(A)]/N[male,b1,app] *
        ## prob that pp was the father of B2
        fec_sa_quant[male,(app + (b2-b1)) |> clamp(A),qq] * recip_TRO_SY[male,b2] *
        ## lucky litter factor - same cohort HS problem
        ifelse(b2==b1,
               lucky_litter_par,1)
    }
  )
  
  #blur into SYLSYL
  Pr_HSP_Pat_SYLSYL <- autoloop(
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, s2=SEXES, y2=SAMPY, lc2=LENGTH_CLASSES,
    SUMOVER=list(b1=POPY,b2= POPY),{
      # calculate ages based on sample and birth year
      a1 <- y1 - b1
      a2 <- y2 - b2
      
      # blur the BB probabilities into SYLSYL
      Pr_HSP_Pat_BB[b1,b2] *
        Pr_A_SYL[a1 |> clamp(A), s1, y1, lc1] * 
        Pr_A_SYL[a2 |> clamp(A), s2, y2, lc2]
    }
  )
  
  ## Grandparent-Grandchild with unseen mother
  #idealized probability calculation
  Pr_GPP_Mat_SYLB <- autoloop(
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, b2=POPY, 
    SUMOVER=list(app=A), {
      ## two prob comp
      bpp <- b2 - app
      (bpp <= PopYear1) *
        # prob that app is age A given b2
        fec_sa[female,app] * recip_TRO_SY[female,b2] *
        # POP prob for GP to be app's parent
        Pr_POP_SYLB[s1,y1,lc1,bpp |> clamp(POPY)]
    }
  )
  
  # blurring probability to sum over potential ages based on length
  Pr_GPP_Mat_SYLSYL1 <- autoloop(
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, s2=SEXES, y2=SAMPY, lc2=LENGTH_CLASSES,
    SUMOVER=list(a2 = A), 
    {
      #calculating gk birth year
      b2 <- y2 - a2
      (b2 >= PopYear1) *
        Pr_GPP_Mat_SYLB[s1, y1, lc1, b2 |> clamp(POPY)] *
        Pr_A_SYL[a2, s2, y2, lc2]
    })
  
  Pr_GPP_Mat_SYLSYL2 <- autoloop(
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, s2=SEXES, y2=SAMPY, lc2=LENGTH_CLASSES, {
      Pr_GPP_Mat_SYLSYL1[s2,y2,lc2,s1,y1,lc1]
    }
  )
  
  ## Grandparent-Grandchild with unseen father
  #idealized probability calculation
  Pr_GPP_Pat_SYLB <- autoloop(
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, b2=POPY, 
    SUMOVER=list(app=A), {
      ## two prob comp
      bpp <- b2 - app
      (bpp <= PopYear1) *
        # prob that app is age A given b2
        fec_sa[male,app] * recip_TRO_SY[male,b2] *
        # POP prob for GP to be app's parent
        Pr_POP_SYLB[s1,y1,lc1,bpp |> clamp(POPY)]
    }
  )
  
  # blurring probability to sum over potential ages based on length
  Pr_GPP_Pat_SYLSYL1 <- autoloop(
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, s2=SEXES, y2=SAMPY, lc2=LENGTH_CLASSES,
    SUMOVER=list(a2= A), 
    {
      #calculating gk birth year
      b2 <- y2 - a2
      (b2 >= PopYear1) *
        Pr_GPP_Pat_SYLB[s1, y1, lc1, b2 |> clamp(POPY)] *
        Pr_A_SYL[a2, s2, y2, lc2]
    })
  
  #for weird case where the grandparent was indiv 2
  Pr_GPP_Pat_SYLSYL2 <- autoloop(
    s1=SEXES, y1=SAMPY, lc1=LENGTH_CLASSES, s2=SEXES, y2=SAMPY, lc2=LENGTH_CLASSES, {
      Pr_GPP_Pat_SYLSYL1[s2,y2,lc2,s1,y1,lc1]
    }
  )
  
  ## combining half-sib and GPP probs
  Pr_TKP_SYLSYL <- Pr_HSP_Mat_SYLSYL + Pr_HSP_Pat_SYLSYL + 
    Pr_GPP_Mat_SYLSYL1 + Pr_GPP_Mat_SYLSYL2 + 
    Pr_GPP_Pat_SYLSYL1 + Pr_GPP_Pat_SYLSYL2
  
  ## generating exp kin pairs
  E_N_TKP_SYLSYL <- n_comps_POP_SYLSYL* Pr_TKP_SYLSYL
  
  ## summing likelihoods for all cases
  nll <- nll - sum(dpois(c(N_TKP_SYLSYL),c(E_N_TKP_SYLSYL),log = T),na.rm = T)
  
  ####
  ### Add mitoDNA for half-sibling groups
  ####
  for (i in 1:nrow(TKPairs)) {
    #grabbing values fron the mtDNA data input
    TKP_s1 <- TKPairs$s1[i]
    TKP_y1 <- TKPairs$y1[i]
    TKP_l1 <- TKPairs$l1[i]
    TKP_s2 <- TKPairs$s2[i]
    TKP_y2 <- TKPairs$y2[i]
    TKP_l2 <- TKPairs$l2[i]
    
    TKP_h1 <- TKPairs$h1[i]
    TKP_h2 <- TKPairs$h2[i]
    
    ## conditionalish probabilities for all four cases
    #HSP maternal - if they match it's 1, no match 0
    Pr_h2_h1_HSP_Mat <- ifelse(TKP_h1==TKP_h2,1,0)
    #HSP paternal - prob of match is just the haplotype frequency
    Pr_h2_h1_HSP_Pat <- HapFreq[TKP_h2]
    
    #GPP maternal - nested ifelse where sex of grandparent is male/female
    #for maternal grandmother, it is the 1/0 prob like mat HS
    # for maternal grandfather, is the hap frequency
    Pr_h2_h1_GPP_Mat1 <- ifelse(TKP_s1==female,ifelse(TKP_h1==TKP_h2,1,0),HapFreq[TKP_h2])
    # case where bc of weird length stuff indiv 2 is the maternal grandparent instead of grandchild
    Pr_h2_h1_GPP_Mat2 <- ifelse(TKP_s2==female,ifelse(TKP_h2==TKP_h1,1,0),HapFreq[TKP_h2])
    #GPP paternal - heritance is broken so whatevs
    Pr_h2_h1_GPP_Pat1 <- HapFreq[TKP_h2]
    Pr_h2_h1_GPP_Pat2 <- HapFreq[TKP_h2]
    ## stitch together 6 cases into the overall probability of h1 and h2 to be added to the likelihood
    Pr_h2_h1_TKP_SYLSYL <- 
      # prob of mat HSP for SYLSYL pair i multiplied by the haplotype probability
      Pr_HSP_Mat_SYLSYL[SLICE=TKP_s1,SLICE=TKP_y1,SLICE=TKP_l1,SLICE=TKP_s2,SLICE=TKP_y2,SLICE=TKP_l2] * Pr_h2_h1_HSP_Mat +
      # prob of pat HSP for SYLSYL pair i multiplied by the haplotype probability
      Pr_HSP_Pat_SYLSYL[SLICE=TKP_s1,SLICE=TKP_y1,SLICE=TKP_l1,SLICE=TKP_s2,SLICE=TKP_y2,SLICE=TKP_l2] * Pr_h2_h1_HSP_Pat +
      # prob of mat GPP for SYLSYL pair i multiplied by the haplotype probability
      Pr_GPP_Mat_SYLSYL1[SLICE=TKP_s1,SLICE=TKP_y1,SLICE=TKP_l1,SLICE=TKP_s2,SLICE=TKP_y2,SLICE=TKP_l2] * Pr_h2_h1_GPP_Mat1 +
      # prob of mat GPP for SYLSYL pair i multiplied by the haplotype probability
      # for the problem case when grandchild is indiv 1
      Pr_GPP_Mat_SYLSYL2[SLICE=TKP_s1,SLICE=TKP_y1,SLICE=TKP_l1,SLICE=TKP_s2,SLICE=TKP_y2,SLICE=TKP_l2] * Pr_h2_h1_GPP_Mat2 +
      # prob of pat GPP for SYLSYL pair i multiplied by the haplotype probability
      Pr_GPP_Pat_SYLSYL1[SLICE=TKP_s1,SLICE=TKP_y1,SLICE=TKP_l1,SLICE=TKP_s2,SLICE=TKP_y2,SLICE=TKP_l2] * Pr_h2_h1_GPP_Pat1 +
      # prob of pat GPP for SYLSYL pair i multiplied by the haplotype probability
      # for the problem case when grandchild is indiv 1
      Pr_GPP_Pat_SYLSYL2[SLICE=TKP_s1,SLICE=TKP_y1,SLICE=TKP_l1,SLICE=TKP_s2,SLICE=TKP_y2,SLICE=TKP_l2] * Pr_h2_h1_GPP_Pat2
    
    #divide probability the same denominator - Pr_TKP_SYLSYL because we've established that they are
    #SOKPs before we started the loop
    full_hap_Pr <- Pr_h2_h1_TKP_SYLSYL/Pr_TKP_SYLSYL[SLICE=TKP_s1,SLICE=TKP_y1,
                                                     SLICE=TKP_l1,SLICE=TKP_s2,
                                                     SLICE=TKP_y2,SLICE=TKP_l2]
    
    #add prob to the likelihood
    nll <- nll - log(full_hap_Pr)
  }
  
  ## Age composition data
  #samp_asyl is an input created for the data
  #adds to the likelihood probabilities based on the age composition
  #list of ages identified for each length bin in a given sample year
  #only use when we have across-sample age data, which we don't have much of rn
  #nll <- nll - sum(log(Pr_A_SYL[MATSUB=samp_asyl]))
  
  REPORTO(N, Z, E_POP_SYLSYL, E_N_TKP_SYLSYL, fec_sa)
  ##Return nll
  nll
})

##FOR NOW DISABLE recruitment parameters
tmpout <- new_f(parm)

tmbmap = list(log_rec=as.factor(rep(NA,length(parm$log_rec))))

testo = MakeADFun(new_f,parm,random=c("rw_log_rec"))

testo$gr()

repp <- testo$report()
opt <- nlminb(testo$par, testo$fn, testo$gr)

sdr <- sdreport(testo)

pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")

