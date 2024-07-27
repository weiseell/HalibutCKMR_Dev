#include <TMB.hpp>

enum sex_t{
  male = 0,
  female = 1
};

template<class Type>
matrix<Type> prob_len_at_age(vector<Type> mean_len_at_age,vector<Type> sd_len_at_age,vector<Type> length_bins,vector<Type> ages){
  int A = ages.size();
  int L = length_bins.size();

  matrix<Type> pla(L,A);

  for(int a = 0; a < A; ++a){
    Type mu_a = mean_len_at_age(a);
    Type sd_a = sd_len_at_age(a);


    for(int l = 0; l < L; ++l){
      if(l == 0){
	pla(l,a) = pnorm(length_bins(l),mu_a,sd_a);
      }else if(l == L-1){
	pla(l,a) = 1-pnorm(length_bins(l-1),mu_a,sd_a);
      }else{
	pla(l,a) = pnorm(length_bins(l),mu_a,sd_a)-pnorm(length_bins(l-1),mu_a,sd_a);
      }
    }


  }


  return pla;
}

template<class Type>
Type make_fecundity(Type alpha, Type bexp, Type len){
  return alpha*pow(len,bexp);
}
  

template<class Type>
Type objective_function<Type>::operator() ()
{

  //Load in Data

  
  //Total number of ages
  DATA_INTEGER(A);
  //Total number of years
  DATA_INTEGER(Y);
  //Total number of length bins
  DATA_INTEGER(L);

  //Columns for sex
  DATA_MATRIX(mean_len_at_age);
  DATA_MATRIX(sd_len_at_age);
  DATA_VECTOR(length_bins);
  DATA_VECTOR(ages);
  
  //The actual observations
  //This is just the stuff from the nonzero data.frame but adjusted so the numbers are correct
  DATA_IVECTOR(sex_x);
  DATA_IVECTOR(AgeAtSamp_x);
  DATA_IVECTOR(SampYear_x);
  DATA_IVECTOR(lbin_ind);

  DATA_IVECTOR(sex_y);
  DATA_IVECTOR(AgeAtSamp_y);
  DATA_IVECTOR(SampYear_y);

  DATA_VECTOR(ncomp);
  DATA_VECTOR(POP);

  //Add HSP, GPGC etc. 
  
  //Load in Parameters
  
  PARAMETER_VECTOR(log_rec);
  PARAMETER(log_rec_sd);
  PARAMETER_VECTOR(log_avg_Z);
  PARAMETER(log_init_abundance);
  //columns for each sex
  PARAMETER_MATRIX(log_fec_parm);


  //Transformations of parameters
  Type rec_sd = exp(log_rec_sd);
  Type init_N = exp(log_init_abundance);
  

  
  Type nll = 0;


  array<Type> N(A,Y,2);
  array<Type> Z(A,Y,2);
  array<Type> P_l_at_age(L,A-1,2);

  
  for(int s = 0; s < 2; ++s){
    vector<Type> mean_len_at_age_s = mean_len_at_age.col(s);
    vector<Type> sd_len_at_age_s = sd_len_at_age.col(s);
    matrix<Type> P_l_at_age_s = prob_len_at_age(mean_len_at_age_s,sd_len_at_age_s,length_bins,ages);
    //.col here will grab the entire remaining matrix
    P_l_at_age.col(s) = P_l_at_age_s.array();
  }
  
  //This will probably be what's actually needed 
  array<Type> fecun_L(L,2);
  //But since the simulation has no concept of length...
  array<Type> fecun_A(A,2);

  for(int s = 0; s < 2; ++s){
    Type alpha = exp(log_fec_parm(0,s));
    Type bexp = exp(log_fec_parm(1,s));
    Type fec_denom = 0;
    fecun_A(0,s) = 0;
    for(int a = 1; a < A; ++a){
      Type runfec = 0;
      for(int i = 0; i < L; ++i){
	Type l = length_bins(i);
	Type fec = make_fecundity(alpha,bexp,l);
	fecun_L(i,s) = fec;
	runfec += P_l_at_age(i,a-1,s)*fec;
	fec_denom += fec_denom+runfec;
      }
      fecun_A(a,s) = runfec;
    }
    
  }
		

  //Random walk for recruitment
  for(int y = 1; y < Y; ++y){
    nll -= dnorm(log_rec(y),log_rec(y-1),rec_sd,true);
  }

  //Put in the recruitment
  for(int y = 0; y < Y; ++y){
    N(0,y,male) = exp(log_rec(y))/2;
    N(0,y,female) = exp(log_rec(y))/2;
  }

  //Put in the mortality
  for(int s = 0; s < 2; ++s){
    for(int a = 0; a < A; ++a){
      for(int y = 0; y < Y; ++y){
	Z(a,y,s) = exp(log_avg_Z(s));
      }
    }
  }

  //Initial abundance
  
  //Cumlative survival calculation
  matrix<Type> cumsurv(A-1,2);
  for(int s = 0; s < 2; ++s){
    cumsurv(0,s) = exp(-Z(0,0,s));
    for(int a = 1; a < A-1; ++a){
      cumsurv(a,s) = exp(-Z(a,0,s));
      cumsurv(a,s) *= cumsurv(a-1,s);
    }
  }

  Type sumcs = cumsurv.sum();

  for(int a = 1; a < A; ++a){
    for(int s = 0; s < 2; ++s){
      N(a,0,s) = init_N*(cumsurv(a-1,s)/sumcs);
    }
  }

  //Do rest of pop dynamics
  for(int y = 1; y < Y; ++y){
    for(int a = 1; a < A; ++a){
      for(int s = 0; s < 2; ++s){
	N(a,y,s) = N(a-1,y-1,s)*exp(-Z(a-1,y-1,s));
      }
    }
  }
  
  //OBSERVATION PART

  //POPs
  int n_pop_obs = POP.size();
  //Some stuff for debuggins
  vector<Type> probs(n_pop_obs);
  vector<Type> nlls(n_pop_obs);
  
  for(int i = 0; i < n_pop_obs; ++i){
    //Grab all the stuff we need
    //Offspring birth year
    int B2 = SampYear_y[i]-AgeAtSamp_y[i];
    //Parent birth year
    int B1 = SampYear_x[i]-AgeAtSamp_x[i];
    //parent sex
    int sex1 = sex_x[i];
    //offspring sex (not really needed but eh)
    int sex2 = sex_y[i];

    Type prob;
    //I guess this is lethal only sampling...
    if(SampYear_x[i] < B2){
      prob = 0;
    }else if(B1 >= B2 || B2 < 0){
      prob = 0;
    }else{
      //PARENTS AGE AT BIRTH, subtract one for zero indexing
      
      int p_a_at_birth = B2-B1-1;
      prob = fecun_A(p_a_at_birth,sex1)/N(p_a_at_birth,B2,sex1);
    }
    probs(i) = prob;
    Type cur_comp = ncomp[i];
    Type cur_pop = POP[i];

    //Now we can just do the nll part
    nlls(i) = dbinom(cur_pop,cur_comp,prob,true);
    nll -= dbinom(cur_pop,cur_comp,prob,true);
  }
  

  
  //HSPs


  //GPGC?


  //Mitochondrial? This is spelt wrong. Do Americans use spelt or spelled?

  //REPORT SOME STUFF
  REPORT(fecun_A);
  REPORT(N);
  REPORT(Z);
  REPORT(cumsurv);
  REPORT(probs);
  REPORT(nlls);
  REPORT(P_l_at_age);
  REPORT(fecun_L);

  return nll;
}
