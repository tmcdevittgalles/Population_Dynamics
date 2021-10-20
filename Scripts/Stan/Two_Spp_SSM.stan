/// Trying to include a second species

data {
  int< lower = 1 > N; // total sample size 
  int< lower = 0 > P; // number of parameters
  int< lower = 0 > Time; // total of number of days
  int <lower = 1 > nSpp; // number of species 
  int Y[N]; // count data
  int Julian[N]; // Julian date of samples
  vector[N] offset; // vector of offsets
  matrix[Time,P] X ; // abiotic covariates
  int State_Julian[Time] ;// Tracking julian date of predictor variables
  matrix[N, nSpp] Spp; // matirx of species identity
}

parameters {
  vector[P] beta_sp1; // parameters 
  vector[P] beta_sp2; // parameters 
  real< lower = 0, upper = 1 > rho;// autoregressive parameter for  each plot
  real< lower = 0 > sigma;
  matrix[Time,nSpp] m;  // predicted states
}

transformed parameters{
  matrix[Time, nSpp] r ;
  
  for(i in 1:Time-1){
    r[i,1] = (X[i,]*beta_sp1);
    r[i,2] = (X[i,]*beta_sp2); 
   }
  
}

model {
  
  // process model
  for( s in 1:nSpp){
     m[1,s] ~ normal(-10, sigma)  ;
     
    for(t in 2:Time) {  // state space model

          (m[t,s]) ~ normal( ((rho)*m[t-1,s]) + r[t-1,s] , sigma);
    }
  }
 
  // data model 
 
    for (i in 1:N) {
      int date;
      date = Julian[i];
      if( Spp[i,1] == 1){
        Y[i] ~ poisson_log( exp(m[date,1]) + (offset[i]));
      }else{
         Y[i] ~ poisson_log( exp(m[date,2]) + (offset[i]));
      }
    }

  // priors
  
    beta_sp1 ~ cauchy(0, 2) ;
    beta_sp2 ~ cauchy(0, 2) ;
    rho ~ beta(8,2);
    sigma ~ cauchy(0, 2) ;
  

}

