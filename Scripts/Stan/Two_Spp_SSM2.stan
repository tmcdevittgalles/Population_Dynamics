/// Trying to include a second species

data {
  int< lower = 1 > N; // total sample size 
  int< lower = 0 > P; // number of parameters
  int< lower = 0 > Time; // total of number of days
  int< lower = 1 > nSpp; // number of species 
  int Y[N]; // count data
  int Julian[N]; // Julian date of samples
  vector[N] offset; // vector of offsets
  matrix[Time*nSpp,P] X ; // abiotic covariates
  matrix[N, nSpp] Spp; // matirx of species identity
}

parameters {
  vector[P] beta; // parameters  
  real< lower = 0, upper = 1 > rho[nSpp];// autoregressive parameter for  each plot
  real< lower = 0 > sigma;
  matrix[Time,nSpp] m;  // predicted states
}

transformed parameters{
  matrix[Time-1, nSpp] r ;
  
  for( s in 1:nSpp){
     int index;
     int index_r;

    if( s == 1 ){
    index = 1;
    }else{
    index = Time +1 ;
    }
    index_r = 1;
    
  for(i in index:(Time * s)-1){
  
    r[index_r,s] = (X[i,]*beta); 
    
    index_r = index_r+1;
   }
  
  }
 
}

model {
  
  // process model
  for( s in 1:nSpp){
     m[1,s] ~ normal(-10, sigma)  ;
    
     
    for(t in 2:Time) {  // state space model

          (m[t,s]) ~ normal( ((rho[s])*m[t-1,s]) + r[t-1,s] , sigma);
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
  
    beta ~ cauchy(0, 2) ;
    rho ~ beta(4,2);
    sigma ~ cauchy(0,5) ;
  

}

