/// Incorporating sampling date and allowing the model to predict counts for 
/// all days and not just sampling day, lets see

data {
  int< lower = 1 > N; // total sample size 
  int< lower = 0 > P; // number of parameters
  int< lower =0 > dP; // number of data model parameters
  int< lower = 0 > Time; // total of number of day
  int Y[N]; // count data
  int Julian[N]; // Julian date of samples
  vector[N] offset; // vector of offsets
  matrix[N, dP] DX ; // Matrix of data model covariates
  matrix[Time,P] X ; // abiotic covariates
}

parameters {
  vector[P] beta; // parameters 
  vector[dP] dBeta; // data model parameters
  real< lower = 0, upper = 1 > rho;// autoregressive parameter for  each plot
  real< lower = 0 > sigma;
  vector[Time] m;  // predicted states
}

model {
  
  // process model
  
     m[1] ~ normal(-10, sigma)  ;
     
  for(t in 2:Time) {  // state space model

        (m[t]) ~ normal( ((rho)*m[t-1]) + (X[t-1,]*beta), sigma);
  }
 
  // data model 
 
    for (i in 1:N) {
      int date;
      date = Julian[i];
        Y[i] ~ poisson_log( exp(m[date]) + DX[i,]*dBeta+ (offset[i]));
    }

  // priors
  
    beta ~ normal(0, 5) ;
    rho ~ beta(8,2);
    sigma ~ exponential(.01) ;
    dBeta ~ normal(0, 5) ;
  

}

