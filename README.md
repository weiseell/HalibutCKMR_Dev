# HalibutCKMR_Dev

Github repository for model development during Codeathon week July 19-Aug 2
* does not contain simulation or kin identification scripts, just TMB/RTMB model development things


At the beginning of the workshop, the Github has:
- CBUT folder: has TMB version of the POP only version of the Halibut model
  - .cpp with the C++ model - needs to be compiled
  - .R with the R code for concatenating the data inputs and parameters to run the C++ model, and optimize the parameters given the data inputs
- Inputs: has simulation data and necessary inputs for the current model versions
  - comparisons and number of POPs and half-siblings
  - Census size from simulation
  - length/age mean and sd from Armsworthy and Campana 2010
- ModelFunction: has the functions that run within the compiled part of the model
  - prob_la: 
  - quant_la: 

