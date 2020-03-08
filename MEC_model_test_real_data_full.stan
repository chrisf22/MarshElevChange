data {
  int<lower=0> N;
  int<lower=0> NN;
  int<lower=0> NNN;
  vector[N] data_vec;
  vector[N] year_index;
  vector[N] restrict;
  int<lower=0,upper=1> veg[N];
  int<lower=1,upper=NN> site_index[N];
  int<lower=1,upper=NNN> patches_index[N];
  int<lower=1,upper=4> geo_index[N];
}

parameters {
  vector[NN] inter;
  vector[NNN] patch;
  real<lower=0> stdev;
  real<lower=0> inter_stdev;
  real<lower=0> patch_stdev;
  real<lower=0> geo_stdev;
  vector[4] geo;
  real B;
  real C;
  real D;
  real inter_mu;
}

transformed parameters {
  vector[N] mu;
  for(i in 1:N){
    mu[i] = inter[site_index[i]] + (D + geo[geo_index[i]] + patch[patches_index[i]] + B*restrict[i] + C*veg[i])*year_index[i];
  }
}

model {
  stdev ~ uniform(0, 100);
  for(q in 1:4){
  geo[q] ~ normal(0, geo_stdev);
  }
  B ~ normal(0, 10);
  C ~ normal(0, 10);
  D ~ normal(0, 20);
  inter_mu ~ normal(0, 500);
  inter_stdev ~ uniform(0, 100);
  patch_stdev ~ uniform(0, 100);
  geo_stdev ~ uniform(0, 10);
  for(z in 1:NN){
  inter[z] ~ normal(inter_mu, site_stdev);
  }
  for(n in 1:NNN){
  patch[n] ~ normal(0, patch_stdev);
  }
  for(i in 1:N){
  data_vec[i] ~ normal(mu[i], stdev);
  }
}
