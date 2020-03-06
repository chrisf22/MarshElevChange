data {
  int<lower=0> N;
  int<lower=0> NN;
  vector[N] data_vec;
  vector[N] year_index;
  int<lower=1,upper=4> region_index[N];
  int<lower=1,upper=NN> site_index[N];
}

parameters {
  vector[NN] inter;
  vector[NN] slope;
  real<lower=0> stdev;
  real<lower=0> inter_stdev;
  real<lower=0> slope_stdev;
  vector[4] inter_mu;
  real slope_mu;
}

transformed parameters {
  vector[N] mu;
  for(i in 1:N){
    mu[i] = inter_mu[region_index[i]] + inter[site_index[i]] + slope[site_index[i]]*year_index[i];
  }
}

model {
  stdev ~ uniform(0, 1000);
  for(q in 1:4){
  inter_mu[q] ~ normal(0, 100);
  }
  slope_mu ~ normal(0, 100);
  inter_stdev ~ uniform(0, 1000);
  slope_stdev ~ uniform(0, 1000);
  for(z in 1:NN){
    inter[z] ~ normal(0, inter_stdev);
    slope[z] ~ normal(slope_mu, slope_stdev);
  }
  for(i in 1:N){
  data_vec[i] ~ normal(mu[i], stdev);
  }
}
