source("Model/DataSet.R", chdir=T)
source("Model/EMParallel2.R", chdir=T)
library("snow")

# # ===============================================
# # Set parameters
# # ===============================================
SEED = 13299
N = 20
M = 50
p = 0.7; q = 0.3; gamma = 2; alpha = 0.03; beta = 0.01
rate_obs = 0.6
process = "ER"
process_est = "ER"
B = 50000
H = 10

# ===============================================
# Generate synthetic data
# ===============================================
param_val = c(p, q, gamma, alpha, beta)
set.seed(SEED)
data = DataSet (N, M, param_val, rate_obs, process = process, start='dense')
yaxis = sapply(data$net_obs_noise, getSizeGCC, order = 1)
cat('data gcc:', yaxis, '\n')

# Data observed
t_obs = data$t_obs # observation moments
Y = data$net_obs_noise 

# Parameter estimation
init = c(0.5, 0.5, 0.5, 0.5, 0.5)
thr = 1e-1
# B = 50000
D = 10
# H = 10
burnIn = ceiling(0.0*M); phi = B*0.9
burnIn_mcmc = ceiling(0.0*D)
MaxIter = 30
prop_max = 50
particleParallel = T
MHParallel = F
seed_parallel = NULL
cat("Settings: p, q, gamma, alpha, beta = ", c(p, q, gamma, alpha, beta), "\n")
cat("          N, M, B, H, D, process, process_est, rate_obs, SEED = ", c(N, M, B, H, D), process, process_est, rate_obs, SEED, "\n")
cat("          Initialization = ", init, "\n")

cat('==================================\n')
cat('Estimating parameters...\n')
cat('==================================\n')
time = proc.time()
results = EM_parallel (Y, t_obs, B, process = process_est, thr, init = init, D = D, MaxIter = MaxIter,
                       prop_max = prop_max, burnIn, phi, H, particleParallel, MHParallel, burnIn_mcmc, seed_parallel)
time_used = proc.time() - time
cat('\n\n')
cat("estimation(p,q,gamma,alpha,beta): \n", results$estimation, "\n")
cat("Time elapsed:", time_used[3],"\n")

cat("%%%%%%%%%%%% RESULTS %%%%%%%%%%%%:\n")
cat("Time elapsed: ", time_used[3],"\n")
cat("Initialization: ", init, "\n")
cat("Truth: ", c(p, q, gamma, alpha, beta), "\n")
cat("Settings: p, q, gamma, alpha, beta = ", c(p, q, gamma, alpha, beta), "\n")
cat("          N, M, B, H, D, process, process_est, rate_obs, SEED = ", c(N, M, B, H, D), process, process_est, rate_obs, SEED, "\n")
cat("(thr, burnIn, burnIn_mcmc, phi, MaxIter): ", c(thr, burnIn, burnIn_mcmc, phi, MaxIter), "\n")
cat("(particleParallel, clusterSizePF):", c(particleParallel, as.numeric(Sys.getenv("NSLOTS"))), "\n")
cat("MHParallel, clusterSizeMH):", c(MHParallel, as.numeric(Sys.getenv("NSLOTS"))), "\n")
cat("Estimates (p, q, gamma, alpha, beta): ", results$estimation, "\n")
cat("iteration: ", results$iteration, "\n")


