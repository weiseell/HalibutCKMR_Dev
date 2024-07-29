#Conversion of Basic CKMR - Brook Trout - to Simplified Halibut Model
#!# VERY ROUGH MODEL!!!!!!!!
## Written by: Ellie Weise
## Last edited: April 12, 2024

#load libraries
library(RTMB)

#load halibut data
load("RTMB_Input/POPs_HS_RTMB_Sim_10000recruits.rda")
load("SimResultsCompiled.rda")
load("Input/SCALmodeldata2022.rdata")
load("Input/HAlibutModelOutputs_NumbersSurvivalSSB.rdata")
la_SCAL <- read.table("Input/lengthage_lvb_fromSCAL.txt",header = T)

##set up data and parameters to give to RTMB model, must be as named lists
## matching what's in the model
dat = list()
#number of sampling years
dat$Y = 6
#maximum age
dat$A = 50
# maximum length for halibut?
dat$L_max = 250

## Calculations needed for dat inputs ####
## data parameters for mortality
#!# total pop data is from the simulation!
#!# simplified survivorship from an 'average year'
#!# will need to find a way to incorporate all halibut data years?
## Survival
S_males<-matrix(NA,dim(N_males)[1],dim(N_males)[2],dimnames = dimnames(N_males))

for(i in 1:(nrow(N_males)-1)){
  for(j in 1:(ncol(N_males)-1)){
    S_males[i,j] <- N_males[i+1,j+1]/N_males[i,j]
  }
  
}

S_females<-matrix(NA,dim(N_females)[1],dim(N_females)[2],dimnames = dimnames(N_females))

for(i in 1:(nrow(N_females)-1)){
  for(j in 1:(ncol(N_females)-1)){
    S_females[i,j] <- N_females[i+1,j+1]/N_females[i,j]
  }
  
}

#calculated survival
Surv_sim_m <- c(1,S_males[1:28,1],rep_len(S_males[28,1],21))
Surv_sim_f <- c(1,S_females[1:28,1],rep_len(S_females[28,1],21))

#calculated mortality
Z_m <- -log(Surv_sim_m)
Z_f <- -log(Surv_sim_f)

#init surviving fraction for each age class based on survival
#init vector based on first two years
init_frac_m <- c(Surv_sim_m[1],Surv_sim_m[2])
init_frac_f <- c(Surv_sim_m[1],Surv_sim_m[2])
#add rest of ages based on who survives to the next year
for (i in 3:dat$A) {
  init_frac_m[i] <- Surv_sim_m[i]*init_frac_m[i-1]
  init_frac_f[i] <- Surv_sim_f[i]*init_frac_f[i-1]
}

#standardizing init_frac so all the years add up to 1 (proportion of total pop)
init_frac_f <- init_frac_f/sum(init_frac_f)
init_frac_m <- init_frac_m/sum(init_frac_m)

#!#set Z and init_frac as input variables

## data parameters for maturity
#!# these are estimated from the SCAL model
#!# currently using age for maturity not length
#!# need to switch this in more finished model
# variable is ppn maturity sex-separated by age class
## Maturity calculations
ppnmature <- data.frame(matrix(nrow = 50, ncol = 2))
colnames(ppnmature) <- c("male","female")

M50F <- 8.5
M95F <- 11.5

mat_a <- 1./(1.+exp(-log(19.)*(1:dat$A-M50F)/(M95F-M50F)));
ppnmature[,2] <- round(mat_a,3)

M50M <- 5
M95M <- 8

mat_b <- 1./(1.+exp(-log(19.)*(1:dat$A-M50M)/(M95M-M50M)));
ppnmature[,1] <- round(mat_b,3)

## data parameters for fecundity
#adding the fake plus group to
#the length at age matrix
#!# this will need to be adjusted when we add the plus group complexity
la_dat <- matrix(0,nrow = dat$A,ncol = 2)

la_dat[,1] <- c(la_SCAL[,1],rep_len(la_SCAL[30,1],length.out = 20))
la_dat[,2] <- c(la_SCAL[,2],rep_len(la_SCAL[30,2],length.out = 20))

## Adding dat parameters to the list ####
#mortality dats
dat$Z = as.matrix(data.frame(Z_m = Z_m, Z_f = Z_f))
dat$init_frac = as.matrix(data.frame(init_frac_m = init_frac_m, init_frac_f = init_frac_f))

#maturity dats
dat$ppnmature <- as.matrix(ppnmature)

#fecundity data
#!# these are from the SCAL model (not)
#!# these are for fecundity (eggs) from a paper i think????
dat$alpha = 0.019
dat$wexp = 3.624

# these are length-weight from the SCAL model, use it so the SR relationship holds
dat$alpha = BiologicalParameters$wt_a
dat$wexp = BiologicalParameters$wt_b


dat$LA = as.matrix(la_dat)

# SR parameters from SCAL
dat$reca = 0.259
dat$recb = 0.209


#POPs
dat$POPs = rels

##Our initial guess for parameters
parm = list()
parm$log_init_tot_N = log(2000000) # start with something closer to reality, current pop is ~8m 
#parm$log_fecundity = log(rep_len(1.1,50))

## function to go into RTMB ####
f <- function(parm){
  getAll(parm,dat)
  
  # id numeric variables
  male <- 1
  female <- 2
  
  #setting numbers matrix
  logN_f <- matrix(0,nrow=A,ncol=Y)
  logN_m <- matrix(0,nrow=A,ncol=Y)
  meanfec_sa <- matrix(0,nrow = A, ncol = 2)
  TROF <- matrix(0,nrow = Y, ncol = 2)
  TROM <- matrix(0,nrow = Y, ncol = 2)
  # fecundity function (will need to be converted to length!)
  get_fec <- function(age,sex){
    # generating the weight at age
    weight = alpha * LA[age,sex] ^ wexp
    #creating a fecundity 'weight' based on proportion mature and weight at age
    fec = (ppnmature[age,sex] * weight)
    fec
  }
  
  # generating fecundity for each age and sex group
  for (s in 1:2) {
    for (a in 1:A) {
      #!# add length calc from skeleton here
      meanfec_sa[a,s] <- get_fec(a,s)
      
    }
  }
  # browser()
  
  #!# convert to relative rather than scaled fecundity
  #meanfec_sa[,1] <- meanfec_sa[,1]/min(meanfec_sa[,1])
  #meanfec_sa[,2] <- meanfec_sa[,2]/min(meanfec_sa[,2])
  
  #setting initial first age
  #this may need to be fucked with too based on recruitment data?
  #need logN for males and females since init fractions are different
  #and TRO needs to be calculated separately
  for (i in 1:A) {
    logN_f[i,1] <- log_init_tot_N + log(init_frac[i,female])
    logN_m[i,1] <- log_init_tot_N + log(init_frac[i,male])
  }
  
  ## Recruitment calculation
  # used to estimate total reproductive output
  #!# need to adjust this for sex as well
  for (s in 1:2) {
    for (y in 2:Y) {
      #set total reproductive output counter
      #calc additive reproductive output
      #uses fecundity and estimated age abund of the year before
      SSB <- sum(meanfec_sa[,2]*exp(logN_f[,y-1])) / 10^9 # convert to kt units
      
      # stock-recruit function to convert SSB to age 1 recruits
      
      TROF[y,s] <- reca * SSB /(1 + recb * SSB) * 10^6 # convert to # of recruits
      logN_f[1,y] <- log(TROF[y,s]/2)
      logN_m[1,y] <- log(TROF[y,s]/2)
      #adjust for mortality
      for (a in 2:A) {
        #!# N will need to be adj for the plus group
        #!# and how recruitment/TRO vary from year to year
        logN_f[a,y] <- logN_f[a-1,y-1]-Z[a-1,2]
        logN_m[a,y] <- logN_m[a-1,y-1]-Z[a-1,1]
      }
    }
  }
  
  ## CKMR Section
  #loop to calculate the TRO/fec for all the years/ages (sex separated)
  #this is used within the loop to estimate the probability for each set of pairs
  #fecundity loop
  for (s1 in 1:2) {
    for (a1 in 1:A) {
      runfec <- 0
      for (l1 in 1:L_max) {
        #use the 'get_fec' equation
        #need to create the 'probability of age given length' input
        #will be for each length bin
        #!# what is Pr_l_given_as??
        runfec <- runfec + Pr_l_given_as() * get_fec[a1,s1]
      }
      meanfec_sa[s1,a1] <- runfec
    }
  }
  
  #TRO loop
  TROgen <- matrix(0,nrow = Y, ncol = 2)
  for (s1 in 1:2) {
    for (y1 in 1:Y) {
      patRO <- 0
      for (a1 in 1:A) {
        patRO <- patRO + exp(logN_f[a1,y1]) * meanfec[a1,s1]
      }
      TROgen[s1,y1] <- patRO
    }
  }
  
  # creating a set of nested for loops for the possibilities in POPs
  # create an array to store all of the probabilities in
  #!# how to create a six-dimensional array in R?
  #!# I can't find tutorials for more than 3 dimensions
  Pr_MOP_SYLAYA <- array()
  # generating the POP probabilities
  for (s1 in 1:2) { # Sex of Parent
    for (y1 in 1:Y) { # Sampling Years for Parents
      for (l1 in ) { # Length of Parent
        for (a1 in 1:A) { # Age of Parent
          for (y2 in 1:Y) { # Sampling Year for Offspring
            for (a2 in 1:A) { # Age of Offspring
              #age of offspring
              B2 <- y2-a2
              #age of parent
              B1 <- y1-a1
              
              #probability of zero scenarios
              if(y1 < B2){
                Prob <- 0
              }
              else if(B1 >= B2) {
                Prob <- 0
              }
              
              else{
                #why do we need L_inf here?? Don't we already have a length/age key?
                L_inf_A <- l1/(1-exp())
                #calculate the age at juvenile birth year
                a_b2 <- a1 = (y1-B2)
                #calculate length at juvenile birth year
                #how to do this???
                l_b2
                #get probability for the pair
                Prob <- get_fec(s1,l_b2)/TROgen[s1,B2]
              }
              Pr_MOP_SYLAYA[s1,y1,l1,a1,y2,a2] <- Prob
            }
          }
        }
        #!# create another array here, four dimensions, length is the number of sampling years
        #!# need to fill with zeros?
        Pr_A_SYL <- array(A,L_inf,2,Y)
        Prob_denom <- 0
        for (a1 in 1:A) {
          #!# again how to make this 'prob given length object??'
          Prob_num <- Pr_l_given_as[a1,s1]*exp(logN_f[a1,y1])
          Prob_denom <- Prob_denom + Prob_num
          Pr_A_SYL[a1,l1,s1,y1] <- Prob_num
        }
        
        for (a1 in 1:A) {
          # what does this do????
          Pr_A_SYL[a1,l1,s1,y1] = Pr_A_SYL[a1,l1,s1,y1]/Prob_denom
        }
      }
    }
  }
  
  ## Addition of information from HSPs
  # I have NO IDEA how this is working
  # basically I just transcribed this bit from the outline
  Pr_OHSP_sbb <- array(2,Y,Y)
  
  for (sp in 1:2) {
    if(sp == 1){Ntmp <- logN_m}
    if(sp == 2){Ntmp <- logN_f}
    for (B1 in 1:Y) {
      for (B2 in 1:Y) {
        Pr_HSP_running <- 0
        for (app in 1:A) {
          Pr_surv <- Ntmp[app + (B2-B1), B2]/Ntmp[app,B1]
          pp_ERO_B2 <- 0
          denom <- 0
          
          for (q in 1:num_qs) {
            #!# what is num_qs and the quantile thing?????
            # what is length_at_age_q?? I think it is an array but not sure how to make it
            #right now it's written as a data frame that is age x length, would need to be sex-separated
            #or could be a 3D array
            Pr_pp_is_qth_q_1 <- 1/(num_qs*get_fec(length_at_age_q[q,app],sp))
            denom <- denom + Pr_pp_is_qth_q_1
            
            ##ERO for B2
            pp_ERO_B2_if_q <- get_fec(length_at_age_q[a,app+B2-B1],sp)
            #ERO for all potential birth years
            #weighted by prob of length at age?
            pp_ERO_B2 <- pp_ERO_B2 + Pr_pp_is_qth_q_1 * pp_ERO_B2_if_q
          }
          # deviding ERO by Pr - denominator across all potential ages
          pp_ERO_B2 <- pp_ERO_B2/denom
          
          # probability of shared parent??? Not exactly sure how all this comes together
          Pr_pp_is_1s_mum <- meanfec_sa[sp,app]/TROgen[sp,B1]
          
          #putting the HSP probability together
          Pr_HSP_running <- Pr_HSP_running + Ntmp[app,B1] * Pr_pp_is_1s_mum * Pr_surv * (pp_ERO_B2/TROgen[sp,B2])
          
          ## compensate for same-cohort HSs (lucky litter)
          #need to make a vector of some kind called lucky_litter_factor
          Pr_HSP_running <- Pr_HSP_running * lucky_litter_factor(sp)
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
  for (s in 1:2) {
    for (y in 1:Y) {
      for (l in 1:length(Lvec)) {
        for (a in 1:A) {
          # neg log-likelihood adjustment for oto ages
          nll = nll - dbinom(n_otos_at_age[a,l,s,y],
                            #total number of otoliths, adjust CKMR Prob based on that 
                            n_otos[l,s,y], Pr_A_SYL[s,y,l],T)
        }
      }
    }
  }
  
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
            #!# this array we'd know both ages in the pair
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
  
  ## adding likelihood based of POP composition
  # need to encompass all parts of the Prob vector
  # this one only has parent length, may need to add uncertain age?
  # also need to adjust these dimensions for what the non-zero comparisons are
  for (s1 in 1:2) {
    for (y1 in 1:Y) {
      for (l1 in 1:L_max) {
        for (a1 in 1:A) {
          for (y2 in 1:Y) {
            for (a2 in 1:A) {
              #!# need arrays for these as well
              n_comps <- n_comps_POP_sylaya[s1,y1,l1,a1,y2,a2]
              n_POPs <- n_POP_sylaya[s1,y1,l1,a1,y2,a2]
              
              ## add these to the likelihood
              # number of pops and lambda is comps weighted by the probability
              # of the pairs as calculated in the above POP prob loop
              nll <- nll - dpois(n_POPs,
                                 n_comps*Pr_MOP_SYLAYA[s1,y1,l1,a1,y2,a2],T)
            }
          }
        }
      }
    }
  }
  
  
  ## RTMB optimization is meant to maximize N and mean fecundity
  N <- exp(logN_f)
  
  REPORT(N)
  REPORT(meanfec_sa)
  
  ADREPORT(logN_f)
  ADREPORT(log_fecundity)
  
  # returning the log likelihood
  nll
}

##Make the objective function with RTMB
obj = MakeADFun(f,parm)
obj$gr()

##optimize the model to get the most likely parameters given the data
opt = nlminb(obj$par,obj$fn,obj$gr)

##get the data out from the model, make confidence intervals
report = obj$report()
sdr = sdreport(obj)
logN = as.data.frame(summary(sdr)[which(rownames(summary(sdr)) == "logN"),])
logN$year = sort(rep(10:19,3))
logN$age = rep(1:3,10)
logN$lower = exp(logN[,1] - 1.96*logN[,2])
logN$upper = exp(logN[,1] + 1.96*logN[,2])
logN$estN = exp(logN[,1])

logFec = as.data.frame(summary(sdr)[which(rownames(summary(sdr)) == "log_fecundity"),])
logFec = logFec[1:3,]
logFec$fec = exp(logFec[,1])
logFec$lower = exp(logFec[,1] - 1.96*logFec[,2])
logFec$upper = exp(logFec[,1] + 1.96*logFec[,2])


##compare with true N and fecundity
trueN = stuff$trueN
logN$trueN = as.vector(trueN)
trueFec = stuff$trueFec

library(ggplot2)
ggplot(logN) + 
  geom_line(aes(x=year,y=estN),color="red") + 
  geom_line(aes(x=year,y=trueN),color="black") + 
  geom_ribbon(aes(x=year,ymin=lower,ymax=upper),alpha=0.3) + 
  facet_wrap(~age)

logFec

