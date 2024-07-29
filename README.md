# HalibutCKMR_Dev

Github repository for model development during Codeathon week July 19-Aug 2
* does not contain simulation or kin identification scripts currently, just TMB/RTMB model development things


At the beginning of the workshop, the Github has:
- CBUT folder: has TMB version of the POP only version of the Halibut model
  - .cpp with the C++ model - needs to be compiled
  - .R with the R code for concatenating the data inputs and parameters to run the C++ model, and optimize the parameters given the data inputs
- Inputs: has simulation data and necessary inputs for the current model versions
  - comparisons and number of POPs and half-siblings
  - Census size from simulation
  - length/age mean and sd from Armsworthy and Campana 2010
- ModelFunction: has the functions that run within the compiled part of the model
  - prob_la: generates probability of length at age, sex separated
  - quant_la: generates a given set of quantile length at age, sex separated
  
Here are some proposed topics for what to cover during the codeathon:
1. Revision of Babyn et al paper
2. Intricacies of the HS version of the model
 - Quantiles for length – maybe apply this to the POPs too?
3. Grandparent/Grandchild problem for HS model
4. Incorporating Incomplete Age Data
5. Post-hoc probs for Age and mitochondrial data - related individuals
6. Uncertainty/simulation testing framework
  - 
7. Ne calculation inside CKMR framework
8. Kinference
9. Next steps
