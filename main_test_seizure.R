source("Model/DataSet.R", chdir=T)
source("Model/EMParallel2.R", chdir=T)
source("Model/GetMargLik2.R", chdir=T)
source("Model/BayesFactor.R", chdir=T)
library(snow)
library(R.matlab)

# ===============================================
# Argument parser for command line input
# ===============================================
args = commandArgs(trailingOnly = TRUE)
# Seed, N, M, p, q, gamma, alpha, beta, rate_obs, process, B, H

# ===============================================
# Get parameters from arguement parser
# ===============================================
SEED = as.numeric(args[1]) 
B = as.numeric(args[2])
H = as.numeric(args[3])
data_path = args[4]
start = as.numeric(args[5])
end = as.numeric(args[6])

# data_path = 'Data/BU1_Seizure2_nets_clean_pad100-10_[4-50]Hz_corr_0_lag_scaled.mat'
# start = 410
# end = 450
# ===============================================
# Load data
# ===============================================
cat('Loading data...\n')
data = readMat(data_path)
Y = list() # data observed
t_obs = numeric(0)
id = 1
for (t in start:end){
  t_obs = c(t_obs, data$nets[,,1]$t[t])
  Y[[id]] = data$nets[,,1]$C[,,t]
  id = id + 1
}
rm(data)

# data information
cat('t_obs: ', t_obs, '\n')
yaxis = sapply(Y, getSizeGCC, order = 1)
cat('data gcc:', yaxis, '\n')
rate_obs = length(t_obs)/(t_obs[length(t_obs)] - t_obs[1])
cat('Observation rate is ', rate_obs, '\n')

M = length(t_obs)
N = dim(Y[[1]])[1]
cat('N = ', N, ' M = ', M, '\n')

# Parameter estimation
#init = c(0.5, 0.5, rate_obs, 0.1, 0.1)
#init = c(0.9, 0.1, rate_obs, 0.01, 0.01)
#init = c(0.9, 0.1, 30, 0.001, 0.01)
#init = c(0.9, 0.1, rate_obs, 0.001, 0.01)
init = c(0.9, 0.1, 10, 0.01, 0.01)
thr = 1e-1
# B = 50000
D = 10
# H = 2
burnIn = ceiling(0.0*M); phi = B*0.9
burnIn_mcmc = ceiling(0.0*D)
MaxIter = 30
prop_max = 50
particleParallel = T
MHParallel_ER = F
MHParallel_PR = T
seed_parallel = SEED

cat("          N, M, B, H, D, rate_obs, seed_parallel = ", c(N, M, B, H, D), rate_obs, SEED, "\n")
cat("                                  Initialization = ", init, "\n")

## under ER
cat('==================================\n')
cat('Estimating parameters under ER...\n')
cat('==================================\n')
time = proc.time()
results = EM_parallel (Y, t_obs, B, process = "ER", thr, init = init, D = D, MaxIter = MaxIter,
                         prop_max = prop_max, burnIn, phi, H, particleParallel, MHParallel_ER, burnIn_mcmc, seed_parallel)
time_used = proc.time() - time
para_est_ER = results$estimation
cat('\n\n')
cat("estimation(p,q,gamma,alpha,beta): \n", results$estimation, "\n")
cat("Time elapsed:", time_used[3],"\n")

## under PR
cat('==================================\n')
cat('Estimating parameters under PR...\n')
cat('==================================\n')
init = para_est_ER
time = proc.time()
results = EM_parallel (Y, t_obs, B, process = "PR", thr, init = init, D = D, MaxIter = MaxIter,
                         prop_max = prop_max, burnIn, phi, H, particleParallel, MHParallel_PR, burnIn_mcmc, seed_parallel)
time_used = proc.time() - time
para_est_PR = results$estimation
cat('\n\n')
cat("estimation(p,q,gamma,alpha,beta): \n", results$estimation, "\n")
cat("Time elapsed:", time_used[3],"\n")

# ===============================================
# Approximation for BF
# ===============================================
cat('para_est_ER=', para_est_ER, '\n')
cat('para_est_PR=', para_est_PR, '\n')

margParallel = T
bf = getBayesFactor(N, Y, t_obs, para_est_ER, para_est_PR, B, margParallel, seed_parallel)
cat("Bayes Factor is: ", bf, '\n')
cat('The truth is', process, '\n')

