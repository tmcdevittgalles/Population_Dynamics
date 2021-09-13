/// Incorporating sampling date and allowing the model to predict counts for 
/// all days and not just sampling day, lets see

data {
  int< lower = 1 > N; // total sample size 
  int< lower = 0 > P; // number of parameters
  int< lower = 0 > Time; // total of number of days
  //int G; // number of group
  int Y[N]; // count data
  int Julian[N]; // Julian date of samples
  vector[N] offset; // vector of offsets
  matrix[Time,P] X ; // abiotic covariates
 // int group_samp[G] ; // sample size for each group
}

parameters {
  vector[P] beta; // parameters 
  real< lower = 0, upper = 1 > rho;// autoregressive parameter for  each plot
  real< lower = 0 > sigma[Time];
 
}

transformed parameters{
  // process model
   vector[Time] m;  // predicted states
   
     m[1]= 1  ;
  for(t in 2:Time) {  // state space model

        m[t] =  (rho*m[t-1]) + (X[t-1,]*beta);

  }
}

model {
  
  
  // data model 
 
    for (i in 1:N) {
      int date;
      date = Julian[i];
        Y[i] ~ poisson_log(m[date] + offset[i]);
    }

  // priors
  
    beta ~ cauchy (0, 5) ;
    rho ~  uniform(0,1); //beta(8,2);
    sigma ~ cauchy (0, 5) ;

}

