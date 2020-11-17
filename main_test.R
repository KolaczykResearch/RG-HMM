source("Model/DataSet.R", chdir=T)
source("Model/EMParallel2.R", chdir=T)
source("Model/GetMargLik2.R", chdir=T)
source("Model/BayesFactor.R", chdir=T)
library("snow")

# ===============================================
# Argument parser for command line input
# ===============================================
args = commandArgs(trailingOnly = TRUE)
# Seed, N, M, p, q, gamma, alpha, beta, rate_obs, process, B, H

# ===============================================
# Get parameters from arguement parser
# ===============================================
SEED = as.numeric(args[1]) 
N = as.numeric(args[2])
M = as.numeric(args[3])
p = as.numeric(args[4]); q = as.numeric(args[5])
gamma = as.numeric(args[6]); alpha = as.numeric(args[7]); beta = as.numeric(args[8])
rate_obs = as.numeric(args[9])
process = args[10]
B = as.numeric(args[11])
H = as.numeric(args[12])

# ===============================================
# Generate synthetic data
# ===============================================
param_val = c(p, q, gamma, alpha, beta)
set.seed(SEED)
data = DataSet (N, M, param_val, rate_obs, process = process)
yaxis = sapply(data$net_obs_noise, getSizeGCC, order = 1)
cat('data gcc:', yaxis, '\n')

# Data observed
t_obs = data$t_obs # observation moments
Y = data$net_obs_noise 

# Parameter estimation
init = c(0.9, 0.1, rate_obs, 0.2, 0.2)
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
seed_parallel = NULL
cat("Settings: p, q, gamma, alpha, beta = ", c(p, q, gamma, alpha, beta), "\n")
cat("          N, M, B, H, D, process, rate_obs, SEED = ", c(N, M, B, H, D), process, rate_obs, SEED, "\n")
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
# init = para_est_ER
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

