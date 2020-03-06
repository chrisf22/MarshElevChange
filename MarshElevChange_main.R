# generate simulated data to test stan model  

# for 100 locations and 10 years
simdims <- mat.or.vec(100, 10)

# draw baseline ME rates for each location
ints <- rnorm(length(simdims[,1]), 0, 0.2)

# draw values for rates of change in ME
slopes <- rnorm(length(simdims[,1]), 3, 0.2)

# combine intercepts and slopes to get 
regeq <- function(x){ints + slopes*x}
site_means <- matrix(unlist(lapply(1:10, FUN = regeq)), nrow=100, ncol=10, byrow=FALSE)

# draw values for variation around rates of change
# why doesn't MARGIN=c(1, 2) work as an alternative to adding drawn values of AV to site_means? 
AV <- apply(site_means, MARGIN=2, FUN=rnorm, mean=0, sd=0.1)
# add valyes for annual variation to mean regression equation
site_means_wvar <- AV + site_means 

# add missingness at random 
# index for which years and locations are missing; prob is the proportion of missing observations
missing <- apply(site_means, MARGIN=2, FUN=rbinom, size=1, prob=0.3)
site_means_wvar[missing==1] <- NA

miss_index <- function(x){which(!is.na(x))}
# get a vector for the year of observation for each observation that is not NA
year_list <- apply(site_means_wvar, MARGIN=1, FUN=miss_index)
year_index <- unlist(year_list)

# get a vector for the site index for each observation that is not NA
num_obs <- lapply(year_list, FUN=length)
num_obs_fun <- function(x){rep(x, num_obs[[x]])}
site_index <- unlist(lapply(1:100, FUN=num_obs_fun))

test_index <- which(!is.na(site_means_wvar), arr.ind=TRUE)
test_index <- test_index[order(test_index[,1]), ]
data_vec <- site_means_wvar[test_index]


# load real set data
setwd("/Users/chrisfield/Dropbox/USFWScontract/accretion_data/final_data/")
MACWA_sets <- read.csv(file = "MACWA_23Feb2020.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, quote="")
DENERR_sets <- read.csv(file = "DENERR_23Feb2020.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, quote="")
EBF_sets <- read.csv(file = "EBF_23Feb2020.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, quote="")
MERI_sets <- read.csv(file = "MERI_2March2020.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, quote="")

region <- rep(1, length(MACWA_sets$pin_ht))
MACWA_sets <- cbind(MACWA_sets, region)

region <- rep(2, length(DENERR_sets$pin_ht))
DENERR_sets <- cbind(DENERR_sets, region)

region <- rep(3, length(EBF_sets$pin_ht))
EBF_sets <- cbind(EBF_sets, region)

region <- rep(4, length(MERI_sets$pin_ht))
MERI_sets <- cbind(MERI_sets, region)

MACWA_sets <- rbind(MACWA_sets, DENERR_sets, EBF_sets, MERI_sets)

MACWA_sets$year <- MACWA_sets$year - min(MACWA_sets$year)

sites_num <- unique(MACWA_sets$site)
for(i in 1:length(sites_num)){
  MACWA_sets$site[MACWA_sets$site==sites_num[i]] <- i
}

data_vec <- as.numeric(MACWA_sets$pin_ht)
year_index <- as.numeric(MACWA_sets$year)
site_index <- as.numeric(MACWA_sets$site)
region_index <- as.numeric(MACWA_sets$region)
N <- length(data_vec)
NN <- max(site_index)

library(rstan)
rstan_options(auto_write = TRUE)
setwd("/Users/chrisfield/Documents/folders/GitProjects/MarshElevChange/")
#setwd("/Users/chrisfield/Desktop/")
input_data <- list(data_vec = data_vec, year_index = year_index, site_index=site_index, N=N, NN=NN, region_index=region_index)
fit_cp <- stan(file='MEC_model_test_real_data.stan', data=input_data, iter=20000, warmup=10000, chains=1, seed=483892929, refresh=1200, control = list(adapt_delta = 0.99))
fit_cp <- stan(file='MEC_model_test_real_data.stan', data=input_data, iter=20000, warmup=10000, chains=1, seed=400, refresh=1200, control = list(adapt_delta = 0.99))
print(fit_cp)
params_cp <- as.data.frame(extract(fit_cp, permuted=FALSE))
readRDS("test_normal.rds")

hist(params_cp$`chain:1.ainter_mu`)
hist(params_cp$`chain:1.aslope_mu`)


hist(params_cp$`chain:1.aslope_mu`, main=" ", xlab="mm/year")
slope_means <- mat.or.vec(max(site_index), 1)
for(i in 1:max(site_index)){
  slope_means[i] <- mean(params_cp[,paste("chain:1.slope[", i, "]", sep="")])
  points(slope_means[i], 0, pch=16, col=rgb(0, 0, 0, 0.4))
}





