source("DataSet.R")
source("GetMargLik.R")

# ===============================================
# Set parameters
# ===============================================
set.seed(10221211)
M = 50
N = 3
p = 0.7; q = 0.3; gamma = 2; alpha = 0.03; beta = 0.01
rate_obs = 0.6
process = "ER"

# ===============================================
# Generate synthetic data
# ===============================================
param_val = c(p, q, gamma, alpha, beta)
data = DataSet (N, M, param_val, rate_obs, process = process)

# Data observed
t_obs = data$t_obs # observation moments
net_obs_noise = data$net_obs_noise 

# para_est = results$estimation
para_est = c(0.7, 0.3, 2, 0.03, 0.01)
nChain = 2000
lenChain = NULL
process = "ER"
prob = getMargLik(N, net_obs_noise, t_obs, nChain, lenChain, para_est, process)
cat("\n", prob)
