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

modelfile <- "2_full.txt"
x <- cbind(intercept,age,male,bmi,chol,
           male*bmi,tx*dose,
           tx*dose*month,tx*dose*age,tx*dose*male,tx*dose*bmi,
           tx*dose*male*bmi)
n_covariates <- dim(x)[2]
tau_b <- 1



data <- list( "vite", "n", "n_patients", "n_covariates", "ptid", "x", "tau_b" )

inits <- function() {
  list( beta = rnorm( n_covariates, 0, 1 ),
        gamma = rnorm( n_patients, 0, 1 ),
        tau_bc = runif( 1, 0, 2 ),
        tau_g = runif( 1, 0, 2 ) )
}

parameters <- c("mu", "beta", "gamma", "sigma_bc", "sigma_g" )

start_time <- Sys.time()
vite_2_full <- bugs( data, inits, parameters,
                          working.directory=bc_path, model.file=modelfile, OpenBUGS.pgm=BUGS_path,
                          n.chains=chains, n.burnin=burn_in, n.iter=iterations+burn_in, n.thin=thin,
                          DIC=TRUE, clearWD = FALSE, debug=FALSE)
end_time <- Sys.time()
end_time - start_time

#
# DIC
vite_2_full$DIC

#
# Diagnostic plots

library(coda)
vite_full <- read.bugs(bc_chainfile)
cumuplot(vite_full[,1:6],probs=c(0.05,0.5,0.95))
cumuplot(vite_full[,7:12],probs=c(0.05,0.5,0.95))

autocorr.plot(vite_full[,1:6], lag.max=100)
autocorr.plot(vite_full[,7:12], lag.max=100)

# Proportions of betas larger than 0
sum(unlist(vite_full[,1]) > 0) / iterations
sum(unlist(vite_full[,2]) > 0) / iterations
sum(unlist(vite_full[,3]) > 0) / iterations
sum(unlist(vite_full[,4]) > 0) / iterations
sum(unlist(vite_full[,5]) > 0) / iterations
sum(unlist(vite_full[,6]) > 0) / iterations
sum(unlist(vite_full[,7]) > 0) / iterations
sum(unlist(vite_full[,8]) > 0) / iterations
sum(unlist(vite_full[,9]) > 0) / iterations
sum(unlist(vite_full[,10]) > 0) / iterations
sum(unlist(vite_full[,11]) > 0) / iterations
sum(unlist(vite_full[,12]) > 0) / iterations

plot(vite_full[,1:3])
plot(vite_full[,4:6])
plot(vite_full[,7:9])
plot(vite_full[,10:12])

plot(vite_full[,737])  # sigma vite
plot(vite_full[,738])  # sigma g

#
# Residual diagnostics

vite_fitted <- c()
for(i in 1:n){
  vite_fitted <- c(vite_fitted, median(vite_2_full$sims.list$mu[,i]))
}

resids <- vite - vite_fitted

plot(density(resids))

shapiro.test(resids)

plot(vite_fitted, resids)



####
# REDUCED MODEL
x <- cbind(intercept,bmi,chol,
           tx*dose,
           tx*dose*month,tx*dose*bmi,
           tx*dose*male*bmi)
n_covariates <- dim(x)[2]


start_time <- Sys.time()
vite_2_red <- bugs( data, inits, parameters,
                         working.directory=bc_path, model.file=modelfile, OpenBUGS.pgm=BUGS_path,
                         n.chains=chains, n.burnin=burn_in, n.iter=iterations+burn_in, n.thin=thin,
                         DIC=TRUE, clearWD = FALSE, debug=FALSE)
end_time <- Sys.time()
end_time - start_time

vite_2_red$DIC

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

vite_red <- read.bugs(bc_chainfile)
cumuplot(vite_red[,1:4],probs=c(0.05,0.5,0.95))
cumuplot(vite_red[,5:7],probs=c(0.05,0.5,0.95))

autocorr.plot(vite_red[,1:4], lag.max=100)
autocorr.plot(vite_red[,5:7], lag.max=100)


# Proportions of betas larger than 0
sum(unlist(vite_red[,1]) > 0) / iterations
sum(unlist(vite_red[,2]) > 0) / iterations
sum(unlist(vite_red[,3]) > 0) / iterations
sum(unlist(vite_red[,4]) > 0) / iterations
sum(unlist(vite_red[,5]) > 0) / iterations
sum(unlist(vite_red[,6]) > 0) / iterations
sum(unlist(vite_red[,7]) > 0) / iterations

plot(vite_red[,1:3])
plot(vite_red[,4:6])
plot(vite_red[,7:9])

## RESIDUALS

vite_fitted_red <- c()
for(i in 1:n){
  vite_fitted_red <- c(vite_fitted_red, median(vite_2_red$sims.list$mu[,i]))
}

resids_vite_red <- vite - vite_fitted_red

plot(density(resids_vite_red))

shapiro.test(resids_vite_red)

plot(vite_fitted_red, resids_vite_red)

# WRITE THE RESIDUALS

#
v <- cbind(month, dose, resids_vite_red, vite_fitted_red, vite)
write.csv(v, "v_resids.csv")


## PROBABILITY INTERVALS AND POINT ESTIMATES

vite_red_mode <- c()
for(i in 1:n_covariates){
  vite_red_mode <- c(vite_red_mode, getmode(vite_2_red$sims.list$beta[,i]))
}

vite_red_median <- c()
for(i in 1:n_covariates){
  vite_red_median <- c(vite_red_median, median(vite_2_red$sims.list$beta[,i]))
}

vite_red_mean <- c()
for(i in 1:n_covariates){
  vite_red_mean <- c(vite_red_mean, mean(vite_2_red$sims.list$beta[,i]))
}

vite_red_mode
vite_red_median
vite_red_mean


# intervals
quantile(vite_2_red$sims.list$beta[,1],probs=c(0.15,0.85))
quantile(vite_2_red$sims.list$beta[,2],probs=c(0.15,0.85))
quantile(vite_2_red$sims.list$beta[,3],probs=c(0.15,0.85))
quantile(vite_2_red$sims.list$beta[,4],probs=c(0.15,0.85))
quantile(vite_2_red$sims.list$beta[,5],probs=c(0.15,0.85))
quantile(vite_2_red$sims.list$beta[,6],probs=c(0.15,0.85))
quantile(vite_2_red$sims.list$beta[,7],probs=c(0.15,0.85))

