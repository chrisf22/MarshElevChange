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