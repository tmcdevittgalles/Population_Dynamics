/**
 * Time Series of Count Data:
 * Notes:
 *   - intercept only (no predictors),
 *   - stationarity requires rho < 1, but not enforced by constraint
 *   - priors introduced here
 */
data {
  int<lower=0> N; // number of time points
  //int<lower = 0 > P; // number of parameters ( not including alpha)
  int y[N]; // Count data
  vector[N] offset; // vector of time of trsts
  
}
parameters {
  //real alpha;           // intercept
  real<lower=0> rho;    // autoregression parameter
  real<lower=0> sigma;  // noise scale on latent state evolution
  vector<lower=0>[N] u;  // Predicted states
  //vector[P] beta; // parameters 
}

model {
  // priors
  //alpha ~ normal(0, 10);
  rho ~ lognormal(0, 2);
  sigma ~ lognormal(0, 5);
  //beta ~ normal(0,10);

  // latent state
  u[1] ~ lognormal(0,2) ;
  for (t in 2:N)
  u[t] ~ normal( rho * (u[t - 1]), sigma);

  // likelihood
  for (t in 1:N)
    y[t] ~ poisson_log( u[t] + log(offset[t]) );
}