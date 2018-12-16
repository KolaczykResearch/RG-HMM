source("DataSet.R")
source("EMParallel.R")

# ===============================================
# Set parameters
# ===============================================
set.seed(102212)
M = 100
N = 50
p = 0.7; q = 0.3; gamma = 10; alpha = 0.03; beta = 0.01
rate_obs = 5

# ===============================================
# Generate synthetic data
# ===============================================
param_val = c(p, q, gamma, alpha, beta)
data = DataSet (N, M, param_val, rate_obs, process = "ER")

# Data observed
t_obs = data$t_obs # observation moments
net_obs_noise = data$net_obs_noise # observed networks with noise
 
# Ground truth in hidden layer
net_obs = data$net_obs # true networks at observational moments 
bin_obs = data$bin_obs # binary variables at observational moments
net_hidden = data$net_hidden # true networks at all transitioning moments 
bin_hidden = data$bin_hidden # binary variables at all transitioning moments
t_all = data$t_all # all transitioning time points 

# ===============================================
# Run EM on synthetic data
# ===============================================
set.seed(1022)
Y = net_obs_noise
# init = c(0.7, 0.3, 2, 0.3, 0.2)
init = c(p, q, gamma, alpha, beta)
B = 50000
thr = 1e-1
D = 100
H = 100
burnIn = ceiling(0.1*M); phi = B*0.9
MaxIter = 30
prop_max = 50
particleParallel = F; # clusterSizePF = 3
MHParallel = T; # clusterSizeMH = 3
# outFile = "output.txt"

cat("Settings: p, q, gamma, alpha, beta = ", c(p, q, gamma, alpha, beta), "\n")
cat("          N, M, B, H, D = ", c(N, M, B, H, D), "\n")
cat("          rate_obs = ", rate_obs, "\n")

time = proc.time()
results = EM_parallel (Y, t_obs, B, process = "ER", thr, init = init, D = D, MaxIter = MaxIter, 
                       prop_max = prop_max, burnIn, phi, H, particleParallel, MHParallel)
time_used = proc.time() - time
# time_used

cat("%%%%%%%%%%%% RESULTS %%%%%%%%%%%%:\n")
cat("Time elapsed:", time_used[3],"\n")
cat("Initialization:", init, "\n")
cat("Truth:", c(p, q, gamma, alpha, beta), "\n")
cat("Settings: (N, M, B, H, D)", c(N, M, B, H, D), "\n")
cat("(thr, burnIn, phi, MaxIter):", c(thr, burnIn, phi, MaxIter), "\n")
cat("(particleParallel, clusterSizePF):", c(particleParallel, as.numeric(Sys.getenv("NSLOTS"))), "\n")
cat("MHParallel, clusterSizeMH):", c(MHParallel, as.numeric(Sys.getenv("NSLOTS"))), "\n")
cat("estimation(p,q,gamma,alpha,beta): \n", results$estimation, "\n")
cat("iteration:\n", results$iteration, "\n")

