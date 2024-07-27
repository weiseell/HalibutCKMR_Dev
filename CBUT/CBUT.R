## load libraries
library(TMB)
library(tidyverse)

#compile C++ section of model
compile("Desktop/HalibutCKMR_Dev/CBUT/CBUT.cpp")
dyn.load(dynlib("Desktop/HalibutCKMR_Dev/CBUT/CBUT"))

#load simulated data inputs
load("Desktop/HalibutCKMR_Dev/Inputs/CensusSize_Sim_071224.rda")
load("Desktop/HalibutCKMR_Dev/Inputs/POPs_df_072224.rda")
la_key <- read.table("Desktop/HalibutCKMR_Dev/Inputs/LA_MeanSD.txt",header = T, sep = "\t",stringsAsFactors = F)

#source functions for RTMB model
source("Desktop/HalibutCKMR_Dev/ModelFunction/prob_la.R")

##Should they be centered around a midpoint?
len_bins = seq(10,220,10)

datt = list()
datt$A = 30
datt$Y = 15
datt$L = length(len_bins)

##I might have sexes mixed up...
f_L_means = la_key |>
    filter(sex == 2) |>
    filter(AgeClass <= 30) |>
    select(mean)

f_L_sds = la_key |>
    filter(sex == 2) |>
    filter(AgeClass <= 30) |>
    select(sd)

m_L_means = la_key |>
    filter(sex == 1) |>
    filter(AgeClass <= 30) |>
    select(mean)

m_L_sds = la_key |>
    filter(sex == 1) |>
    filter(AgeClass <= 30) |>
    select(sd)


datt$mean_len_at_age = matrix(c(m_L_means[,1],f_L_means[,1]),ncol=2)
datt$sd_len_at_age = matrix(c(m_L_sds[,1],f_L_sds[,1]),ncol=2)
datt$length_bins = len_bins
##This is just for the PLA
datt$ages = 2:datt$A

###We can't use values where the juvenile is born more than 15 years in the past because it's outside what the model considers...
nonzero = nonzero |>
    mutate(j_by = SampYear.y - AgeAtSamp.y) |>
    filter(j_by >= 0)

##Adjust for zero indexing
datt$sex_x = nonzero$sex.x-1
datt$AgeAtSamp_x = nonzero$AgeAtSamp.x-1
datt$SampYear_x = nonzero$SampYear.x-1
datt$lbin_ind = nonzero$lbin-1

datt$sex_y = nonzero$sex.y-1
datt$AgeAtSamp_y = nonzero$AgeAtSamp.y-1
datt$SampYear_y = nonzero$SampYear.y-1

datt$ncomp = nonzero$ncomp
datt$POP = nonzero$POP

Narray = array(CS_181_195$abundance,c(50,15,2))

parm = list()
parm$log_rec = log(Narray[1,,1] + Narray[1,,2])
parm$log_rec_sd = log(0.3)
parm$log_init_abundance <- log(sum(Narray[-1,1,]))

parm$log_avg_Z = log(c(0.1,0.1))
parm$log_fec_parm = log(matrix(c(3.624,0.019,3.624,0.019),nrow=2))

##Turn off recruitement
mapp = list()
mapp$log_rec = as.factor(rep(NA,length(parm$log_rec)))
mapp$log_rec_sd = as.factor(NA)

system.time(obj <- MakeADFun(datt,parm,mapp,random="log_rec"))
repp = obj$report()

opt = nlminb(obj$par,obj$fn,obj$gr)
