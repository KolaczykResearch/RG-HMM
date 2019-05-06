# """
# This module generates synthetic network data and "snapshot" observations  
# """

source("utilities.R")
source("GetSimSeq.R")

# DataSet()--------------------------------------
# --------------------Input: --------------------
# N: Number of nodes.
# M: Total number of observational time points.
# param_val: true values of p, q, gamma, alpha, beta
# rate_obs: sampling rate for observations, 
#           eg. rate_obs = 2: 2 obs. per unit time, i.e 1 obs. every 1/2 unit time
# process: "ER" or "PR"
# --------------------Output: -------------------
# -------------------
# Synthetic data ###
# -------------------
# net_obs_noise: Network sequence at observational moments with noise. 
# t_obs: Observational time points.
# ---------------------
# Ground truth part ###
# ---------------------
# net_obs: Network sequence at observational moments with no noise, i.e. 
#              A list of length(t_obs) elements, each containing N tims N adjacency matrix.
# net_hidden: Network sequence at all transitioning moments, including net_start
# bin_obs: Binary variable sequence at observational moments. 
# bin_hidden: Binary variable sequence at all transitioning moments, including w0 = 1 by default.
# t_all: All transitioning moments, including 0 for net_start. 
DataSet = function(N, M, param_val, rate_obs, process = "ER"){
  p = param_val[1]
  q = param_val[2]
  gamma = param_val[3]
  alpha = param_val[4]
  beta = param_val[5]
  t_obs = cumsum(rep(1/rate_obs, M))
  # initialize the starting network with a random graph with 1/5 of edges
  net_start = matrix(0, N, N)
  edges = sample_modified(get_triangular(1:N^2, N), ceiling(choose(N,2)/100))
  for (i in 1:length(edges)){
    pos = get_pos(edges[i], N)
    net_start[pos[1], pos[2]] = 1
    net_start[pos[2], pos[1]] = 1
  }
  # generate sequence of networks from process "ER" or "PR"
  res = get_SimSeq_all(p, q, gamma, N, t_obs, net_start = net_start, process = process)
  net_obs = res$network_seq_obs
  net_obs_noise = add_Noise(net_obs, alpha, beta)
  return(list(net_obs_noise = net_obs_noise, t_obs = t_obs, 
              net_obs = net_obs, bin_obs = res$bin_var_obs,
              net_hidden = res$network_seq_hidden, bin_hidden = res$bin_var_hidden,
              t_all = res$t_transition))
}

# ===========================================
# ===========================================
# N = 5; M = 10
# p = 0.7; q = 0.3; gamma = 2; alpha = 0.3; beta = 0.2
# rate_obs = 0.5
# param_val = c(p, q, gamma, alpha, beta)
# data = DataSet (N, M, param_val, rate_obs, process = "ER")
# 
# data$t_obs
# data$net_obs_noise
# 
# data$net_obs
# data$bin_obs
# data$bin_hidden

