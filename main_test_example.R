source("Model/DataSet.R", chdir=T)
source("Model/EMParallel2.R", chdir=T)
source("Model/GetMargLik2.R", chdir=T)
source("Model/BayesFactor.R", chdir=T)
library("snow")

# # ===============================================
# # Set parameters
# # ===============================================
SEED = 13299
N = 10
M = 14
p = 0.9; q = 0.1; gamma = 2; alpha = 0.01; beta = 0.01
rate_obs = 1.5
process = "ER"
B = 50000
H = 4
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
# H = 10
burnIn = ceiling(0.0*M); phi = B*0.9
burnIn_mcmc = ceiling(0.0*D)
MaxIter = 30
prop_max = 50
particleParallel = T
MHParallel_ER = F
MHParallel_PR = T
seed_parallel = 13299
cat("Settings: p, q, gamma, alpha, beta = ", c(p, q, gamma, alpha, beta), "\n")
cat("          N, M, B, H, D, process, rate_obs, SEED = ", c(N, M, B, H, D), process, rate_obs, SEED, "\n")

## under ER
cat('==================================\n')
cat('Estimating parameters under ER...\n')
cat('==================================\n')
time = proc.time()
results = EM_parallel (Y, t_obs, B, process = "ER", thr, init = init, D = D, MaxIter = MaxIter,
                       prop_max = prop_max, burnIn, phi, H, particleParallel, MHParallel_ER, burnIn_mcmc)
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
                       prop_max = prop_max, burnIn, phi, H, particleParallel, MHParallel_PR, burnIn_mcmc)
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


