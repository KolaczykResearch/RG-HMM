source("DataSet.R", chdir = T)
source("MHSamples.R")
source("EMParallel.R")

test_mcmc = function(Y, bin_obs, t_obs, process = "ER", D = 50, prop_max = 100,
                     MHParallel = FALSE, burnIn_mcmc = 0, burnIn = 0, H = 1, p_cur, q_cur, gamma_cur){
  cat("--> Given a sequence of true networks, start E-step...\n")
  time1 = proc.time()
  # Sample H true network in the hidden layer from B particles
  # then for each of them sample D paths that connects the truth
  # \hat p = E[E[\sum(I{w_{n-1} = 0, w_m = 1)]]/E[E[\sum(I_{w_{n-1} = 0})]]
  counts_vec = rep(0, 5)
  M = length(t_obs)
  ids = 1
  partic = list()
  partic$W = matrix(bin_obs, 1, M)
  particles = array(list(), 1)
  particles[[1]] = net_obs
  partic$particles = particles

  if(MHParallel){
    cat("   Sampling MH paths for b = ...\n   ")
    clusterSize = as.numeric(Sys.getenv("NSLOTS"))
    if (is.na(clusterSize)) clusterSize = 1
    #initialize a cluster for parallel computing
    cluster = makeCluster(clusterSize, type = "SOCK")
    #sample D paths that connects Y_true, W_true
    results = clusterApply(cluster, 1:(H*(M-1-burnIn)), getSamplePath_parallel, D, partic, t_obs,
                           p_cur, q_cur, gamma_cur, process = process, prop_max, ids, burnIn, burnIn_mcmc)
    # counts_vec = c(sum_01, sum_0, sum_10, sum_1, sum_R)
    for (cl in 1:(H*(M-1-burnIn))){
      counts_vec = counts_vec + results[[cl]]
    }
    stopCluster(cluster)
  }else{
    cat("   Sampling MH paths for m = ...\n   ")
    for (h in 1:H){
      cat(h, ":", ids[h], "...,")
      b = ids[h]
      Y_true = partic$particles[[b]]
      W_true = partic$W[b, ]
      #sample D paths that connects Y_true, W_true
      for (m in (burnIn+1):(M-1)){
        results = getSamplePath_nonparallel(m, D, Y_true, W_true, t_obs,
                                            p_cur, q_cur, gamma_cur, process = process, prop_max, burnIn_mcmc)
        counts_vec = counts_vec + results
      }
    }
  }
  time_used2 = proc.time() - time1
  cat("  Finished E-step.\n")
  cat("  Time elapsed:", time_used2[3],"\n")
  pqgamma = update_pqgamma(counts_vec, H, D, t_obs, M)
  p_nex = pqgamma[1]
  q_nex = pqgamma[2]
  gamma_nex = pqgamma[3]
  return(list(p_nex = p_nex, q_nex = q_nex, gamma_nex = gamma_nex))
}

# test_mcmc = function(Y, bin_obs, t_obs, process = "ER", D = 50, prop_max = 100, 
#                      MHParallel = FALSE, burnIn_mcmc = 0, burnIn = 0, H = 1, p_cur, q_cur, gamma_cur){
#   cat("--> Given a sequence of true networks, start E-step...\n")
#   time1 = proc.time()
#   # Sample H true network in the hidden layer from B particles
#   # then for each of them sample D paths that connects the truth
#   # \hat p = E[E[\sum(I{w_{n-1} = 0, w_m = 1)]]/E[E[\sum(I_{w_{n-1} = 0})]]
#   counts_vec = rep(0, 5)
#   M = length(t_obs)
#   ids = 1
#   partic = list()
#   partic$W = matrix(bin_obs, 1, M)
#   particles = array(list(), 1)
#   particles[[1]] = net_obs
#   partic$particles = particles
#   
#   cat("   Sampling MH paths for b = ...\n   ")
#   for (h in 1:H){
#     cat(h, ":", ids[h], "...,")
#     b = ids[h]
#     Y_true = partic$particles[[b]]
#     W_true = partic$W[b,]
#     # sample D paths that connects the truth, (parallel for M here, not D)
#     if (MHParallel){
#       clusterSize = as.numeric(Sys.getenv("NSLOTS"))
#       if (is.na(clusterSize)) clusterSize = 1
#       #initialize a cluster for parallel computing
#       cluster = makeCluster(clusterSize, type = "SOCK")
#       #sample D paths that connects Y_true, W_true
#       results = clusterApply(cluster, 1:(M-1), getSamplePath_parallel, D, Y_true, W_true, t_obs,
#                              p_cur, q_cur, gamma_cur, process = process, prop_max)
#       # counts_vec = c(sum_01, sum_0, sum_10, sum_1, sum_R)
#       for (cl in 1:(M-1)){
#         counts_vec = counts_vec + results[[cl]]
#       }
#       stopCluster(cluster)
#     }
#     else{
#       #sample D paths that connects Y_true, W_true
#       for (m in 1:(M-1)){
#         results = getSamplePath_parallel(m, D, Y_true, W_true, t_obs,
#                                          p_cur, q_cur, gamma_cur, process = process, prop_max)
#         counts_vec = counts_vec + results
#       }
#     }
#   }
# 
#   time_used2 = proc.time() - time1
#   cat("  Finished E-step.\n")
#   cat("  Time elapsed:", time_used2[3],"\n")
#   pqgamma = update_pqgamma(counts_vec, H, D, t_obs, M)
#   p_nex = pqgamma[1]
#   q_nex = pqgamma[2]
#   gamma_nex = pqgamma[3]
#   return(list(p_nex = p_nex, q_nex = q_nex, gamma_nex = gamma_nex))
# }



# ===============================================
# Set parameters
# ===============================================
set.seed(10221211)
M = 50
N = 30
p = 0.7; q = 0.3; gamma = 2; alpha = 0.03; beta = 0.02
rate_obs = 0.6
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


# =====================================================
# Estimate p, q, gamma from the sequence of 
# truth networks and binary variables: net_obs, bin_obs
# ====================================================
p_cur = p; q_cur = q; gamma_cur = gamma
test_mcmc(Y, bin_obs, t_obs, process = "ER", D = 10, prop_max = 30, 
          MHParallel = FALSE, burnIn_mcmc = 1, burnIn = 0, H = 1, p_cur, q_cur, gamma_cur)









# 
# 
# 
# 
# ss = getSamplePath(200, y1, w1, y2, w2, t1, t2, p = 0.3, q = 0.7, gamma = 10)
# ss$rate_accept
# 
# t1 = t_obs[3]
# t2 = t_obs[4]
# y1 = net_obs[[3]]
# y2 = net_obs[[4]]
# w1 = bin_obs[3]
# w2 = bin_obs[4]
# t_trans = t_all[t_all > t1 & t_all <= t2]
# y_trans = net_hidden[t_all > t1 & t_all <= t2]
# w_trans = bin_hidden[t_all > t1 & t_all <= t2]
# 
# y1
# y2
# w1; w2
# t1; t2
# t_trans
# y_trans
# w_trans
# length(t_trans)
# # make visualization
# # vis_input = function(){}
# 
# 
# netSeq_path = get_shortestPath(y1, w1, y2, w2)$netSeq_path
# binSeq_path = get_shortestPath(y1, w1, y2, w2)$binSeq_path
# w1
# w2 
# 




