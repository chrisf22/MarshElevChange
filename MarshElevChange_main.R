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


# load set data
setwd("/Users/chrisfield/Dropbox/USFWScontract/accretion_data/final_data/")
MACWA_sets <- read.csv(file = "MACWA_23Feb2020.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, quote="")
DENERR_sets <- read.csv(file = "DENERR_23Feb2020.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, quote="")
EBF_sets <- read.csv(file = "EBF_23Feb2020.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, quote="")
MERI_sets <- read.csv(file = "MERI_2March2020.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, quote="")
sets_site_list <- read.csv(file = "set_sites_all_SHARP_geo1_wveg_forR.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, quote="")

#create a vector to index which dataset the rows originate from (to keep track after they are merged)
# MACWA
region <- rep(1, length(MACWA_sets$pin_ht))
MACWA_sets <- cbind(MACWA_sets, region)
# DENERR
region <- rep(2, length(DENERR_sets$pin_ht))
DENERR_sets <- cbind(DENERR_sets, region)
# EBF
region <- rep(3, length(EBF_sets$pin_ht))
EBF_sets <- cbind(EBF_sets, region)
# MERI
region <- rep(4, length(MERI_sets$pin_ht))
MERI_sets <- cbind(MERI_sets, region)

# merge the SET datasets; keep the name for MACWA to match the rest of the code
MACWA_sets <- rbind(MACWA_sets, DENERR_sets, EBF_sets, MERI_sets)

# subtract the minimum year from the year records so that each is the number of years since the first SET measurement
#MACWA_sets$year <- MACWA_sets$year - min(MACWA_sets$year)

# get a list of the unique site names
sites_num <- unique(MACWA_sets$site)

###THIS IS PROBABLY WHERE YOU GET RID OF THE FIRST YEARS###
# go through and standardize years to a patch specific first date, then get rid of all of the rows with the lowest date; loop through site num
starting_year <- rep(0, length(MACWA_sets$year))
MACWA_sets <- cbind(MACWA_sets, starting_year)
for(i in 1:length(sites_num)){
  temp <- MACWA_sets[MACWA_sets$site == sites_num[i], ]
  first_year <- min(temp$year)
  baseline <- temp$pin_ht[temp$year==min(temp$year)]
  MACWA_sets[MACWA_sets$site == sites_num[i], ]$pin_ht <- MACWA_sets[MACWA_sets$site == sites_num[i], ]$pin_ht - baseline
  MACWA_sets[MACWA_sets$site == sites_num[i], ]$year <- MACWA_sets[MACWA_sets$site == sites_num[i], ]$year - first_year
  MACWA_sets[MACWA_sets$site == sites_num[i], ]$starting_year <- first_year
}
MACWA_sets <- MACWA_sets[MACWA_sets$year != 0, ]

# filter the list of SET data to keep only sites that are in one of the datasets
sets_site_list <- sets_site_list[unlist(lapply(sets_site_list$site, function(r) any(r %in% sites_num))), ]

# find any sites in the datasets that are not present in the site list 
# remove rows for the site that is missing from the site list from the data (for now)
MACWA_sets <- MACWA_sets[MACWA_sets$site!=setdiff(sites_num, sets_site_list$site), ]
# for both the datasets and the list of SET data, remove rows for the SETs that are outside of SHARP patches (for now)
MACWA_sets <- MACWA_sets[unlist(lapply(MACWA_sets$site, function(r) any(r %in% sets_site_list$site[!is.na(sets_site_list$PatchID)]))), ]
sets_site_list <- sets_site_list[unlist(lapply(sets_site_list$site, function(r) any(r %in% sets_site_list$site[!is.na(sets_site_list$PatchID)]))), ]


# find the unique sites from the site list 
sites_unique <- unique(sets_site_list$site)
# for both the datasets and the list of SET data, make the site indices numeric
for(i in 1:length(sites_unique)){
  MACWA_sets$site[MACWA_sets$site==sites_unique[i]] <- i
  sets_site_list$site[sets_site_list$site==sites_unique[i]] <- i
}

# find the unique patch IDs from the site list
patches_unique <- unique(sets_site_list$PatchID)
# for both the datasets and the list of SET data, make the patch IDs numeric
for(i in 1:length(patches_unique)){
  sets_site_list$PatchID[sets_site_list$PatchID==patches_unique[i]] <- i
}

# find the unique values for geomorphic setting
geo_unique <- unique(sets_site_list$geo1)
# for both the datasets and the list of SET data, make the indices for geomorphic setting numeric
for(i in 1:length(geo_unique)){
  sets_site_list$geo1[sets_site_list$geo1==geo_unique[i]] <- i
}

# create vectors for Stan 
# index for site
site_index <- as.numeric(MACWA_sets$site)
# index for geomorphic setting
geo_index <- mat.or.vec(length(MACWA_sets$pin_ht), 1)
# index for patch ID
patches_index <- mat.or.vec(length(MACWA_sets$pin_ht), 1)
# index for proportion of patch that is restricted
restrict <- mat.or.vec(length(MACWA_sets$pin_ht), 1)
# index for marsh vegetation type (1 = low marsh)
veg <- mat.or.vec(length(MACWA_sets$pin_ht), 1)
# expand vectors so that they are the length of the data vector, instead of the length of the number of sites (NN vs N in Stan)
for(i in 1:length(MACWA_sets$pin_ht)){
  geo_index[i] <- as.numeric(sets_site_list$geo1[site_index[i]])
  patches_index[i] <- as.numeric(sets_site_list$PatchID[site_index[i]])
  restrict[i] <- as.numeric(sets_site_list$Restric_ha[site_index[i]])/(as.numeric(sets_site_list$area_ha[site_index[i]]))
  veg[i] <- as.numeric(sets_site_list$Veg_cover[site_index[i]])
}

# if there is a NA for restrict, assume there is no restriction
restrict[is.na(restrict)] <- 0
# specify marsh categories that are higher elevation than high marsh as high marsh 
veg[veg!=1&veg!=2] <- 1
# recode indext so that high marsh = 0 and low marsh = 1
veg <- veg - 1

# data vector of elevation measurements
data_vec <- as.numeric(MACWA_sets$pin_ht)
# data for the year of the measurement (in number of years since the earliest measurement)
year_index <- as.numeric(MACWA_sets$year)
#region_index <- as.numeric(MACWA_sets$region)
# number of total measurements
N <- length(data_vec)
# number of sites
NN <- max(site_index)
# number of patches
NNN <- max(patches_index)

#length(geo_index)
#length(patches_index)
#length(restrict)
#length(data_vec)
#length(year_index)
#N
#NN

library(rstan)
rstan_options(auto_write = TRUE)
setwd("/Users/chrisfield/Documents/folders/GitProjects/MarshElevChange/")
#setwd("/Users/chrisfield/Desktop/")
#input_data <- list(data_vec = data_vec, year_index = year_index, site_index=site_index, N=N, NN=NN, region_index=region_index)
input_data <- list(data_vec = data_vec, year_index = year_index, N=N, NN=NN, NNN=NNN, geo_index=geo_index, patches_index=patches_index, site_index=site_index, restrict=restrict, veg=veg)
fit_cp <- stan(file='MEC_model_test_real_data_full.stan', data=input_data, iter=20000, warmup=10000, chains=1, seed=483892929, refresh=1200, control = list(adapt_delta = 0.99))
fit_cp <- stan(file='MEC_model_test_real_data.stan', data=input_data, iter=20000, warmup=10000, chains=1, seed=400, refresh=1200, control = list(adapt_delta = 0.99))
print(fit_cp)
params_cp <- as.data.frame(extract(fit_cp, permuted=FALSE))
readRDS("test_normal.rds")

#hist(params_cp$`chain:1.ainter_mu`)
#hist(params_cp$`chain:1.aslope_mu`)

# create a vector for geomorphic setting of length = the number of patches to feed into prediction equation
geo_pred <- mat.or.vec(length(patches_unique), 1)
# create a vector for slr rate of length = the number of patches to feed into prediction equation
slr <- mat.or.vec(length(patches_unique), 1)
# create a vector for proportion of forest loss in the patch of length = the number of patches to feed into prediction equation
forest <- mat.or.vec(length(patches_unique), 1)
# create a matrix to store patch-level predictions, with patches as columns and steps in the MCMC chain as rows
preds <- mat.or.vec(length(params_cp$`chain:1.D`), length(patches_unique))
# same as the matrix above but predictions are relative to the SLR rate
preds_relative <- mat.or.vec(length(params_cp$`chain:1.D`), length(patches_unique))
for(i in 1:length(patches_unique)){
geo_pred[i] <- as.numeric(sets_site_list$geo1[sets_site_list$PatchID==i][1])
slr[i] <- as.numeric(sets_site_list$noaa_slt[sets_site_list$PatchID==i][1])
forest[i] <- as.numeric(sets_site_list$Loss_prop[sets_site_list$PatchID==i][1])
preds[, i] <- params_cp$`chain:1.D` + params_cp[,paste("chain:1.geo[", geo_pred[i], "]", sep="")] + params_cp$`chain:1.B`*0 + params_cp$`chain:1.C`*0 + params_cp[,paste("chain:1.patch[", i, "]", sep="")]
preds_relative[, i] <- preds[, i] - slr[i]
}

# find the posterior distribution for the proportion of patches that are keeping pace with SLR
prop_above0 <- mat.or.vec(length(params_cp$`chain:1.D`), 1)
# for loop to cycle through each step of the MCMC chain
for(i in 1:length(params_cp$`chain:1.D`)){
  # subtract the slr rate for each column (patch) for each step of the MCMC chain
  prop_above0[i] <- length(which((preds[i, ] - slr) > 0))/length(patches_unique)
}

quartz.options(width=6, height=4)
layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))

par(mar=c(3, 3, 1, 1))
plot(slr, colMeans(preds), ylim=c(-8, 10), xlim=c(2, 5), pch=16, col=rgb(0, 0, 0, 0.8), cex.axis=0.75, xlab=" ", ylab=" ")
abline(coef = c(0,1))
for(i in 1:length(slr)){
segments(slr[i], quantile(preds[,i], c(.025)), slr[i], quantile(preds[,i], c(.975)), lend="butt", lwd=10, rgb(0, 0, 0, 0.1))
}
mtext(side=2, line=2, "Marsh ele. change (mm/year)", cex=0.75)
mtext(side=1, line=2, "Sea level rise (mm/year)", cex=0.75)

par(mar=c(3, 3, 1, 1))
hist(prop_above0, cex.axis=0.75, main=" ")
mtext(side=1, line=2, "Prop. of sites above SLR rate", cex=0.75)

par(mar=c(3, 3, 1, 1))
plot(forest, colMeans(preds_relative), ylim=c(-11, 8), pch=16, col=rgb(0, 0, 0, 0.8), cex.axis=0.75)
for(i in 1:length(slr)){
  segments(forest[i], quantile(preds_relative[,i], c(.025)), forest[i], quantile(preds_relative[,i], c(.975)), lend="butt", lwd=10, rgb(0, 0, 0, 0.1))
}
mtext(side=2, line=2, "Marsh ele. change (mm/year)", cex=0.75)
mtext(side=1, line=2, "Forest loss (% of patch)", cex=0.75)


#########EXTRA CODE###########
hist(params_cp$`chain:1.aslope_mu`, main=" ", xlab="mm/year")
slope_means <- mat.or.vec(max(site_index), 1)
for(i in 1:max(site_index)){
  slope_means[i] <- mean(params_cp[,paste("chain:1.slope[", i, "]", sep="")])
  points(slope_means[i], 0, pch=16, col=rgb(0, 0, 0, 0.4))
}





