data {
  vector[688] data_vec;
  vector[688] year_index;
  int<lower=1,upper=100> site_index[688];
}

parameters {
  real B;
  vector[100] inter;
  real<lower=0> stdev;
}

transformed parameters {
  vector[688] mu;
  for(i in 1:688){
    mu[i] = inter[site_index[i]] + B*year_index[i];
  }
}

model {
  stdev ~ uniform(0, 100);
  B ~ normal(0, 100);
  for(z in 1:100){
    inter[z] ~ normal(0, 100);
  }
  for(i in 1:100){
  data_vec[i] ~ normal(mu[i], stdev);
  }
}
