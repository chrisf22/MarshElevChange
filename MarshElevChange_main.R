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
# mute this when using a baseline of zero for each location instead of separate intercepts
#MACWA_sets$year <- MACWA_sets$year - min(MACWA_sets$year)

# get a list of the unique site names
sites_num <- unique(MACWA_sets$site)

###THIS IS PROBABLY WHERE YOU GET RID OF THE FIRST YEARS###
# go through and standardize years to a patch specific first date, then get rid of all of the rows with the lowest date; loop through site num
# mute this for loop and command when using separate intercepts for locations instead of standardizing to zero
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

# export list of sites used in the final analysis for map
sets_site_out <- sets_site_list[, c("site", "lat", "long")]
write.csv(sets_site_out, "/Users/chrisfield/Dropbox/USFWScontract/accretion_data/final_data/SETs_post_analysis.csv")

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
slt <- mat.or.vec(length(MACWA_sets$pin_ht), 1)
# expand vectors so that they are the length of the data vector, instead of the length of the number of sites (NN vs N in Stan)
for(i in 1:length(MACWA_sets$pin_ht)){
  geo_index[i] <- as.numeric(sets_site_list$geo1[site_index[i]])
  patches_index[i] <- as.numeric(sets_site_list$PatchID[site_index[i]])
  restrict[i] <- as.numeric(sets_site_list$Restric_ha[site_index[i]])/(as.numeric(sets_site_list$area_ha[site_index[i]]))
  veg[i] <- as.numeric(sets_site_list$Veg_cover[site_index[i]])
  slt[i] <- as.numeric(sets_site_list$noaa_slt[site_index[i]])
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
input_data <- list(data_vec = data_vec, year_index = year_index, N=N, NN=NN, NNN=NNN, geo_index=geo_index, patches_index=patches_index, site_index=site_index, restrict=restrict, veg=veg, slt=slt)
#input_data <- list(data_vec = data_vec, year_index = year_index, N=N, NNN=NNN, geo_index=geo_index, patches_index=patches_index, restrict=restrict, veg=veg)
fit_cp <- stan(file='MEC_model_test_real_data_full_relative.stan', data=input_data, iter=20000, warmup=10000, chains=1, seed=483892929, refresh=1200, control = list(adapt_delta = 0.99))
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
preds[, i] <- params_cp$`chain:1.D` + params_cp[,paste("chain:1.geo[", geo_pred[i], "]", sep="")] + params_cp$`chain:1.B`*0 + params_cp$`chain:1.C`*0 + params_cp[,paste("chain:1.patch[", i, "]", sep="")] + params_cp$`chain:1.E`*slr[i]
preds_relative[, i] <- preds[, i] - slr[i]
}

slt_x <- seq(1.5, 6, by=0.05)
slt_line <- mat.or.vec(length(params_cp$`chain:1.D`), length(slt_x))
upper <- mat.or.vec(length(slt_x), 1)
lower <- mat.or.vec(length(slt_x), 1)
for(i in 1:length(slt_x)){
slt_line[, i] <- params_cp$`chain:1.D` + params_cp$`chain:1.E`*slt_x[i]
upper[i] <- quantile(slt_line[, i], c(0.975))
lower[i] <- quantile(slt_line[, i], c(0.025))
}


# find the posterior distribution for the proportion of patches that are keeping pace with SLR
prop_above0 <- mat.or.vec(length(params_cp$`chain:1.D`), 1)
forest_corr <- mat.or.vec(length(params_cp$`chain:1.D`), 1)
# for loop to cycle through each step of the MCMC chain
for(i in 1:length(params_cp$`chain:1.D`)){
  # subtract the slr rate for each column (patch) for each step of the MCMC chain
  prop_above0[i] <- length(which((preds[i, ] - slr) > 0))/length(patches_unique)
  forest_corr[i] <- cor(forest, preds[i, ])
}

# multi-panel plot to compare MEC with forest loss and slr
quartz.options(width=5.14, height=4)
layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))

par(mar=c(3, 3, 1, 1))
plot(forest, colMeans(preds_relative), ylim=c(-11, 9), xlim=c(-0.005, 0.055), pch=16, col=rgb(0, 0, 0, 0.8), cex.axis=0.75, mgp=c(3, 0.5, 0))
for(i in 1:length(slr)){
  segments(forest[i], quantile(preds_relative[,i], c(.025)), forest[i], quantile(preds_relative[,i], c(.975)), lend="butt", lwd=10, rgb(0, 0, 0, 0.1))
}
mtext(side=2, line=1.75, "Marsh ele. change (mm/year)", cex=0.65)
mtext(side=1, line=1.5, "Forest loss (prop. of patch)", cex=0.65)
abline(lm(colMeans(preds_relative) ~ forest))

par(mar=c(3, 1, 1, 1))
hist(forest_corr, cex.axis=0.75, main=" ", yaxt="n", mgp=c(3, 0.5, 0))
mtext(side=1, line=1.5, "Correlation between forest loss and MEC", cex=0.65)

par(mar=c(3, 3, 1, 1))
plot(slr, colMeans(preds), ylim=c(-8, 12), xlim=c(2.5, 4.55), pch=16, col=rgb(0, 0, 0, 0.8), cex.axis=0.75, xlab=" ", ylab=" ", mgp=c(3, 0.5, 0))
abline(coef = c(0,1), lwd=2)
for(i in 1:length(slr)){
segments(slr[i], quantile(preds[,i], c(.025)), slr[i], quantile(preds[,i], c(.975)), lend="butt", lwd=10, rgb(0, 0, 0, 0.1))
}
mtext(side=2, line=1.75, "Marsh ele. change (mm/year)", cex=0.65)
mtext(side=1, line=1.5, "Sea level rise (mm/year)", cex=0.65)
#lines(slt_x, lower, lwd=2, col=rgb(0, 0, 0, 0.5))
#lines(slt_x, upper, lwd=2, col=rgb(0, 0, 0, 0.5))
#lines(slt_x, colMeans(slt_line), lwd=2, col=rgb(0, 0, 0, 0.5))
abline(lm(colMeans(preds)~slr), col=rgb(0, 0, 0, 0.5), lwd=2)
conf_int <- predict(lm(colMeans(preds)~slr), newdata=data.frame(slr=slt_x), interval="confidence", level = 0.95)
lines(slt_x, conf_int[,2], col=rgb(0, 0, 0, 0.5), lwd=1)
lines(slt_x, conf_int[,3], col=rgb(0, 0, 0, 0.5), lwd=1)

par(mar=c(3, 1, 1, 1))
hist(prop_above0, cex.axis=0.75, main=" ", yaxt="n", mgp=c(3, 0.5, 0))
mtext(side=1, line=1.5, "Prop. of sites above SLR rate", cex=0.65)


# plot to demonstrate variation at multiple levels in a hierarchical model
# create a matrix of posterior predictions for site-level trends, with a row for each step of the MCMC chain
site_means <- mat.or.vec(length(params_cp$`chain:1.D`), max(site_index))
for(i in 1:max(site_index)){
  site_means[, i] <- mean(params_cp[,paste("chain:1.site[", i, "]", sep="")] + params_cp$`chain:1.D` + params_cp$`chain:1.E`*mean(slt))
}

# create a vector of the mean posterior predictions for patch-level trends
patch_means <- mat.or.vec(max(patches_index), 1)
for(i in 1:max(patches_index)){
  patch_means[i] <- mean(params_cp[,paste("chain:1.patch[", i, "]", sep="")] + params_cp$`chain:1.D` + params_cp$`chain:1.E`*mean(slt))
}

# create a vector of the mean posterior prediction for geo-morphic setting level trends
geo_means <- mat.or.vec(max(geo_index), 1)
for(i in 1:max(geo_index)){
  geo_means[i] <- mean(params_cp[,paste("chain:1.geo[", i, "]", sep="")] + params_cp$`chain:1.D` + params_cp$`chain:1.E`*mean(slt))
}

quartz.options(width=5.14, height=4)
par(mar=c(3, 3, 1, 1))
plot(colMeans(site_means), rep(0, max(site_index)), ylim=c(0, 3.1), pch=16, col=rgb(0, 0, 0, 0), cex.axis=0.75, xlab=" ", ylab=" ", bty="n", yaxt="n", mgp=c(3, 0.5, 0))

# add a vertical line for the mean SLR rate in the region
abline(v=mean(unique(slt)), lwd=2, col="cadetblue")

abline(h=0, col=rgb(0, 0, 0, 0.5))
abline(h=1, col=rgb(0, 0, 0, 0.5))
abline(h=2, col=rgb(0, 0, 0, 0.5))
abline(h=3, col=rgb(0, 0, 0, 0.5))

# choose which Patch ID to use as an example for the figure
ex <- 1

# plot the mean and 95% prediction intervals (using mean estimates of the SD) for the sites within a patch
segments(patch_means[ex], 0, patch_means[ex], 1, col="red", lwd=2)
polygon(x=c((patch_means[ex] - mean(params_cp$`chain:1.site_stdev`)*1.98), (patch_means[ex] - mean(params_cp$`chain:1.site_stdev`)*1.98), (patch_means[ex] + mean(params_cp$`chain:1.site_stdev`)*1.98), (patch_means[ex] + mean(params_cp$`chain:1.site_stdev`)*1.98)),
        y=c(0, 1, 1, 0), border=NA, col=rgb(1, 0, 0, 0.2))

# plot the mean and 95% prediction intervals (using mean estimates of the SD) for the patches within a geomorphic setting
segments(geo_means[as.numeric(sets_site_list$geo1[sets_site_list$PatchID==ex][1])], 1, geo_means[as.numeric(sets_site_list$geo1[sets_site_list$PatchID==ex][1])], 2, col="red", lwd=2)
geoplot <- geo_means[as.numeric(sets_site_list$geo1[sets_site_list$PatchID==ex][1])]
polygon(x=c((geoplot - mean(params_cp$`chain:1.patch_stdev`)*1.98), (geoplot - mean(params_cp$`chain:1.patch_stdev`)*1.98), (geoplot + mean(params_cp$`chain:1.patch_stdev`)*1.98), (geoplot + mean(params_cp$`chain:1.patch_stdev`)*1.98)),
        y=c(1, 2, 2, 1), border=NA, col=rgb(1, 0, 0, 0.2))

# plot the mean and 95% prediction intervals (using mean estimates of the SD) for the variation in geomorphic setting from the overall regional mean
hypermean <- mean(params_cp$`chain:1.D` + params_cp$`chain:1.E`*mean(slt))
segments(hypermean, 2, hypermean, 3, col="red", lwd=2)
polygon(x=c((hypermean - median(params_cp$`chain:1.geo_stdev`)*1.98), (hypermean - median(params_cp$`chain:1.geo_stdev`)*1.98), (hypermean + median(params_cp$`chain:1.geo_stdev`)*1.98), (hypermean + median(params_cp$`chain:1.geo_stdev`)*1.98)),
        y=c(2, 3, 3, 2), border=NA, col=rgb(1, 0, 0, 0.2))

# plot the posterior means for site, patch, and geomorphic setting-level predictions
points(patch_means, rep(1, max(patches_index)), col=rgb(0, 0, 0, 0.8), pch=16)
points(geo_means, rep(2, max(geo_index)), col=rgb(0, 0, 0, 0.8), pch=16)
points(mean(params_cp$`chain:1.D` + params_cp$`chain:1.E`*mean(slt)), 3, col="red", pch=16)
points(colMeans(site_means), rep(0, max(site_index)), col=rgb(0, 0, 0, 0.8), pch=16)

# add red markers to the sites, patches, and geomoprhic setting associated with the example (ex)
points(colMeans(site_means[ ,sets_site_list$PatchID==ex]), rep(0, length(colMeans(site_means[ ,sets_site_list$PatchID==ex]))), col="red", pch=16)
points(patch_means[ex], 1, col="red", pch=16)
points(geo_means[as.numeric(sets_site_list$geo1[sets_site_list$PatchID==ex][1])], 2, col="red", pch=16)

mtext(line=1.5, side=1, "Marsh elevation change (mm/year)", cex=0.75)
text(-7.75, 0.1, "SET locations", pos=4, cex=0.75)
text(-7.75, 1.1, "Marsh patches", pos=4, cex=0.75)
text(-7.75, 2.1, "Geomorphic setting", pos=4, cex=0.75)
text(-7.75, 3.1, "Regional mean", pos=4, cex=0.75)


