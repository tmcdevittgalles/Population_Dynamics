/// Trying to include a second species

data {
  int< lower = 1 > N; // total sample size 
  int< lower = 0 > P; // number of parameters
  int< lower = 0 > Time; // total of number of days
  //int< lower = 1 > nSpp; // number of species 
  int Y[N]; // count data
  int Julian[N]; // Julian date of samples
  vector[N] offset; // vector of offsets
  matrix[Time,P] X ; // abiotic covariates
  //matrix[N, nSpp] Spp; // matirx of species identity
}

parameters {
  vector[P] beta; // parameters  
  real <lower =0, upper = 1> rho;// autoregressive parameter for  each plot
  real< lower = 0 > sigma;
  vector[Time] m;  // predicted states
}

//transformed parameters{
  //vector[Time-1] r ;
  

    
//  for(i in 1:(Time)-1){
  
  //  r[i] = (X[i,]*beta); 
    
   
   //}
  
 // }
 


model {
  
  // process model

     m[1] ~ normal(-10, sigma)  ;
    
     
    for(t in 2:Time) {  // state space model

          (m[t]) ~ normal( ((rho)*m[t-1]) +  (X[t-1,]*beta) , sigma);
    }
  
 
  // data model 
 
    for (i in 1:N) {
      int date;
      date = Julian[i];
     
        Y[i] ~ poisson_log( exp(m[date]) + (offset[i]));

    }

  // priors
  
    beta ~ cauchy(0, 2) ;
    rho ~ beta(4,2);
    sigma ~ cauchy(0,5) ;
  

}

