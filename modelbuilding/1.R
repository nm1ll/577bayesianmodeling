options(scipen=999)
#
######### PREPARATORY CODE ########
#
library("R2OpenBUGS")
bc_path <- "C:\\Users\\nmil\\Dropbox\\_UNM_\\STAT 577 Introduction to Bayesian Modeling\\Final Project\\Code"
BUGS_path <- "C:\\Program Files (x86)\\OpenBUGS\\OpenBUGS323\\OpenBUGS.exe"
bc_datafile <- "bcarotene.csv"
bc_chainfile <- "CODAchain1.txt"
setwd(bc_path)

##  Set MCMC parameters.
chains <- 1           ##  One chain of MCMC values is obtained
burn_in <- 10000       ##  The first 1000 iterates from the MCMC chain are discarded
iterations <- 10000   ##  10,000 values are sampled from the posterior distribution
thin <- 1             ##  No thinning is performed

##  Import beta-carotene data an store it in new variable names for use with OpenBUGS.
bc_data <- read.csv(bc_datafile)
ptid <- bc_data$ptid
month <- bc_data$month
bcarot <- bc_data$bcarot
vite <- bc_data$vite
dose <- bc_data$dose
age <- bc_data$age
male <- bc_data$male
bmi <- bc_data$bmi
chol <- bc_data$chol
cauc <- bc_data$cauc
vauc <- bc_data$vauc

##  Define new variables for OpenBUGS to use.
n <- dim(bc_data)[1]
n_patients <- length( unique(ptid) )
intercept <- rep(1,n)
tx <- intercept - (month<4)


#############################################
#############################################
## 1 FULL ##

modelfile <- "1_full.txt"
x <- cbind(intercept,age,male,bmi,chol,
           male*bmi,tx*dose,
           tx*dose*month,(tx*dose*month)^2,
           tx*dose*age,tx*dose*male,tx*dose*bmi,
           tx*dose*male*bmi)
n_covariates <- dim(x)[2]
tau_b <- 1



data <- list( "bcarot", "n", "n_patients", "n_covariates", "ptid", "x", "tau_b" )

inits <- function() {
  list( beta = rnorm( n_covariates, 0, 1 ),
        gamma = rnorm( n_patients, 0, 1 ),
        tau_bc = runif( 1, 0, 2 ),
        tau_g = runif( 1, 0, 2 ) )
}

parameters <- c("mu", "beta", "gamma", "sigma_bc", "sigma_g" )

start_time <- Sys.time()
bcarotene_1_full <- bugs( data, inits, parameters,
                       working.directory=bc_path, model.file=modelfile, OpenBUGS.pgm=BUGS_path,
                       n.chains=chains, n.burnin=burn_in, n.iter=iterations+burn_in, n.thin=thin,
                       DIC=TRUE, clearWD = FALSE, debug=FALSE)
end_time <- Sys.time()
end_time - start_time

#
# DIC
bcarotene_1_full$DIC

#
# Diagnostic plots

library(coda)
bc_full <- read.bugs(bc_chainfile)
cumuplot(bc_full[,1:6],probs=c(0.05,0.5,0.95))
cumuplot(bc_full[,7:12],probs=c(0.05,0.5,0.95))
cumuplot(bc_full[,13:14],probs=c(0.05,0.5,0.95))

autocorr.plot(bc_full[,1:6], lag.max=100)
autocorr.plot(bc_full[,7:12], lag.max=100)
autocorr.plot(bc_full[,13:14], lag.max=100)

# Proportions of betas larger than 0
sum(unlist(bc_full[,1]) > 0) / iterations
sum(unlist(bc_full[,2]) > 0) / iterations
sum(unlist(bc_full[,3]) > 0) / iterations
sum(unlist(bc_full[,4]) > 0) / iterations
sum(unlist(bc_full[,5]) > 0) / iterations
sum(unlist(bc_full[,6]) > 0) / iterations
sum(unlist(bc_full[,7]) > 0) / iterations
sum(unlist(bc_full[,8]) > 0) / iterations
sum(unlist(bc_full[,9]) > 0) / iterations
sum(unlist(bc_full[,10]) > 0) / iterations
sum(unlist(bc_full[,11]) > 0) / iterations
sum(unlist(bc_full[,12]) > 0) / iterations
sum(unlist(bc_full[,13]) > 0) / iterations
sum(unlist(bc_full[,14]) > 0) / iterations

plot(bc_full[,1:3])
plot(bc_full[,4:6])
plot(bc_full[,7:9])
plot(bc_full[,10:12])
plot(bc_full[,13:14])

plot(bc_full[,737])  # sigma bc
plot(bc_full[,738])  # sigma g

#
# Residual diagnostics

bc_full_fitted <- c()
for(i in 1:n){
  bc_full_fitted <- c(bc_full_fitted, median(bcarotene_1_full$sims.list$mu[,i]))
}

resids_bc_full <- bcarot - bc_full_fitted

par(mfrow=c(1,1))
plot(density(resids_bc_full))

shapiro.test(resids_bc_full)

plot(bc_full_fitted, resids_bc_full)




###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

# REDUCED MODEL
x <- cbind(intercept,male,
           tx*dose,
           tx*dose*month,(tx*dose*month)^2,
           tx*dose*age,
           tx*dose*male*bmi)
n_covariates <- dim(x)[2]
tau_b <- 5

start_time <- Sys.time()
bcarotene_1_red <- bugs( data, inits, parameters,
                          working.directory=bc_path, model.file=modelfile, OpenBUGS.pgm=BUGS_path,
                          n.chains=chains, n.burnin=burn_in, n.iter=iterations+burn_in, n.thin=thin,
                          DIC=TRUE, clearWD = FALSE, debug=FALSE)
end_time <- Sys.time()
end_time - start_time


bcarotene_1_red$DIC

###########
# PLOTS 

bc_red <- read.bugs(bc_chainfile)
cumuplot(bc_red[,1:4],probs=c(0.05,0.5,0.95))
cumuplot(bc_red[,5:7],probs=c(0.05,0.5,0.95))

autocorr.plot(bc_red[,1:4], lag.max=100)
autocorr.plot(bc_red[,5:7], lag.max=100)

# Proportions of betas larger than 0
sum(unlist(bc_red[,1]) > 0) / iterations
sum(unlist(bc_red[,2]) > 0) / iterations
sum(unlist(bc_red[,3]) > 0) / iterations
sum(unlist(bc_red[,4]) > 0) / iterations
sum(unlist(bc_red[,5]) > 0) / iterations
sum(unlist(bc_red[,6]) > 0) / iterations
sum(unlist(bc_red[,7]) > 0) / iterations

plot(bc_red[,1:3])
plot(bc_red[,4:6])
plot(bc_red[,7:9])

plot(bc_red[,737])  # sigma bc
plot(bc_red[,738])  # sigma g

#
# Residual diagnostics

bc_red_fitted <- c()
for(i in 1:n){
  bc_red_fitted <- c(bc_red_fitted, median(bcarotene_1_red$sims.list$mu[,i]))
}

resids_bc_red <- bcarot - bc_red_fitted

par(mfrow=c(1,1))
plot(density(resids_bc_red))

shapiro.test(resids_bc_red)

plot(bc_red_fitted, resids_bc_red)

# WRITE THE RESIDUALS

#
bc <- cbind(month, dose, resids_bc_red, bc_red_fitted,bcarot)
write.csv(bc, "bc_resids.csv")


# Intervals

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

bc_red_mode <- c()
for(i in 1:n_covariates){
  bc_red_mode <- c(bc_red_mode, getmode(bcarotene_1_red$sims.list$beta[,i]))
}

bc_red_median <- c()
for(i in 1:n_covariates){
  bc_red_median <- c(bc_red_median, median(bcarotene_1_red$sims.list$beta[,i]))
}

bc_red_mean <- c()
for(i in 1:n_covariates){
  bc_red_mean <- c(bc_red_mean, mean(bcarotene_1_red$sims.list$beta[,i]))
}

bc_red_mode
bc_red_median
bc_red_mean

# intervals
quantile(bcarotene_1_red$sims.list$beta[,1],probs=c(0.15,0.85))
quantile(bcarotene_1_red$sims.list$beta[,2],probs=c(0.15,0.85))
quantile(bcarotene_1_red$sims.list$beta[,3],probs=c(0.15,0.85))
quantile(bcarotene_1_red$sims.list$beta[,4],probs=c(0.15,0.85))
quantile(bcarotene_1_red$sims.list$beta[,5],probs=c(0.15,0.85))
quantile(bcarotene_1_red$sims.list$beta[,6],probs=c(0.15,0.85))
quantile(bcarotene_1_red$sims.list$beta[,7],probs=c(0.15,0.85))



