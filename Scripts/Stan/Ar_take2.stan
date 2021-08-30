data {
  int< lower = 1 > N; // total sample size 
  int< lower = 0 > P; // number of parameters
  int G; // number of group
  int Y[N]; // count data
  vector[N] offset; // vector of offsets
  matrix[N,P] X ; // abiotic covariates
  int group_samp[G] ; // sample size for each group
}

parameters {
  vector[P] beta; // parameters 
  real< lower = 0, upper = 1 > rho[G];// autoregressive parameter for  each plot
  real< lower = 0 > sigma[G];
  vector< lower = 0 >[N] m;  // predicted states
}

model {
  
  // process model
  
  int pos;
  pos = 1;
  for(i in 1:G) {  // loop over groups
    int local_N;
    vector[group_samp[i]] local_y;
    local_N = group_samp[i];
    m[pos] ~ cauchy (0, 1) ;
    
    for (n in pos+1:local_N+pos-1) { // loop over observations
        m[n] ~ lognormal( (rho[i]*m[n-1]) +(X[n-1,]*beta), sigma[i]);
      }
    pos = pos + local_N;
  }
 
  // data model 
 
    for (i in 1:N) {
        Y[i] ~ poisson_log(m[i] + (offset[i]));
    }

  // priors
  
    beta ~ cauchy (0, 5) ;
    rho ~ uniform(0,1);
    sigma ~ cauchy (0, 5) ;

}

