source("DataSet.R")
source("EMParallel2.R")

# ===============================================
# Set parameters
# ===============================================
set.seed(1210112277)
M = 50
N = 20
p = 0.7; q = 0.3; gamma = 2; alpha = 0.03; beta = 0.01
rate_obs = 0.6
init = c(p, q, gamma, alpha, beta)
#init = c(0.1, 0.1, 1, 0.01, 0.01)
B = 50000
n_rep = 100

one_run = function(M, N, p, q, gamma, alpha, beta, rate_obs, init, B, verbose = F){
  # ===============================================
  # Generate synthetic data
  # ===============================================
  param_val = c(p, q, gamma, alpha, beta)
  data = DataSet (N, M, param_val, rate_obs, process = "ER")
  # Data observed
  t_obs = data$t_obs # observation moments
  Y = data$net_obs_noise # observed networks with noise
  # ===============================================
  # Run EM on synthetic data
  # ===============================================
  thr = 1e-1
  D = 10
  H = 10
  burnIn = ceiling(0.0*M); phi = B*0.9
  burnIn_mcmc = ceiling(0.0*D)
  MaxIter = 30
  prop_max = 50
  particleParallel = F; # clusterSizePF = 3
  MHParallel = F; # clusterSizeMH = 3

  if (verbose){
    cat("========================================================================\n")
    cat("Settings: p, q, gamma, alpha, beta = ", c(p, q, gamma, alpha, beta), "\n")
    cat("          N, M, B, H, D = ", c(N, M, B, H, D), "\n")
    cat("          rate_obs = ", rate_obs, "\n")
    cat("          burnIn, burnIn_mcmc = ", c(burnIn, burnIn_mcmc), "\n")
  }
  
  time = proc.time()
  results = EM_parallel (Y, t_obs, B, process = "ER", thr, init = init, D = D, MaxIter = MaxIter, 
                         prop_max = prop_max, burnIn, phi, H, particleParallel, MHParallel, burnIn_mcmc)
  time_used = proc.time() - time
  # time_used

  cat("%%%%%%%%%%%% RESULTS %%%%%%%%%%%%:\n")
  cat("Time elapsed:", time_used[3],"\n")
  cat("Initialization:", init, "\n")
  cat("Truth:", c(p, q, gamma, alpha, beta), "\n")
  cat("Settings: (N, M, B, H, D)", c(N, M, B, H, D), "\n")
  cat("(thr, burnIn, burnIn_mcmc, phi, MaxIter):", c(thr, burnIn, burnIn_mcmc, phi, MaxIter), "\n")
  cat("(particleParallel, clusterSizePF):", c(particleParallel, as.numeric(Sys.getenv("NSLOTS"))), "\n")
  cat("MHParallel, clusterSizeMH):", c(MHParallel, as.numeric(Sys.getenv("NSLOTS"))), "\n")
  cat("estimation(p,q,gamma,alpha,beta): \n", results$estimation, "\n")
  cat("iteration:\n", results$iteration, "\n")
  return (c(results$estimation, results$iteration))
}
ss = replicate(n_rep, one_run(M, N, p, q, gamma, alpha, beta, rate_obs, init, B, verbose = F))
cat("Mean of parameter estimates: \n p, q, gamma, alpha, beta, iter \n", 
    apply(ss, 1, mean))
cat("Standard deviation  of parameter estimates: \n", apply(ss, 1, sd))




