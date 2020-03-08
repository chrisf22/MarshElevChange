data {
  int<lower=0> N;
  int<lower=0> NNN;
  int<lower=0> NN;
  vector[N] data_vec;
  vector[N] year_index;
  vector[N] restrict;
  vector[N] slt;
  int<lower=0,upper=1> veg[N];
  int<lower=1,upper=NN> site_index[N];
  int<lower=1,upper=NNN> patches_index[N];
  int<lower=1,upper=4> geo_index[N];
}

parameters {
  vector[NN] site;
  vector[NNN] patch;
  real<lower=0> stdev;
  real<lower=0> patch_stdev;
  real<lower=0> geo_stdev;
  real<lower=0> site_stdev;
  vector[4] geo;
  real B;
  real C;
  real D;
  real E;
}

transformed parameters {
  vector[N] mu;
  for(i in 1:N){
    mu[i] = (D + geo[geo_index[i]] + patch[patches_index[i]] + B*restrict[i] + C*veg[i] + E*slt[i])*year_index[i];
  }
}

model {
  stdev ~ uniform(0, 100);
  for(q in 1:4){
  geo[q] ~ normal(0, geo_stdev);
  }
  B ~ normal(0, 100);
  C ~ normal(0, 100);
  D ~ normal(0, 100);
  E ~ normal(0, 100);
  patch_stdev ~ uniform(0, 100);
  site_stdev ~ uniform(0, 100);
  geo_stdev ~ uniform(0, 10);
  for(n in 1:NNN){
  patch[n] ~ normal(0, patch_stdev);
  }
  for(z in 1:NN){
  site[z] ~ normal(0, site_stdev);
  }
  for(i in 1:N){
  data_vec[i] ~ normal(mu[i], stdev);
  }
}
