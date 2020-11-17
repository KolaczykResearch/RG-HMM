source("utilities.R")
source("prob_fun.R")
source("GetSimSeq.R")
source("ParticleFiltering.R")

# getMargLik2()-----------------------
# Calculate marginal log-likelihood of observed network sequence 
# N: number of nodes
# Y: observed sequence
# t_obs: obseration moments
# para_est: parameter estimations
# B: number of particles needed for approximation
# parallel: whether parallel or not in generating B particles at each obs time points
getMargLik2 = function(N, Y, t_obs, para_est, B, process, parallel, seed_parallel=NULL){
  M = length(t_obs)
  W = matrix(0, B, M) # each row is the seq of binary var at t_obs
  W[,1] = rep(1, B) # set binary var to be 1 at the first obs. moments
  marglik_log = log(1)
  edgeIndex = matrix(1:(N^2), N, N)[upper.tri(matrix(1:(N^2), N, N))]

  # create B particles at each t_obs, store sequential particles at t_obs in a list 
  # particles is a list of list, particles[[1]],...,particles[[B]]
  particles = array(list(), B)
  Y1 = reformat(Y[[1]], N)
  for (b in 1:B){
    particles[[b]] = list()
    # particles[[b]][[1]] = Y[[1]]
    particles[[b]][[1]] = Y1
  }
  
  for (m in 2:M){
    cat("m = ", m, "...")
    # check weights
    prob_cond = sapply(particles, FUN = function(list){
      res = getCondProb(observe = Y[[m-1]], true = format_back(list[[m-1]], N), para_est[4], para_est[5])
      return (res)})
    # calculate Z_m-1 
    if (m > 2){
      # average = sum(exp(prob_cond))/B
      average_log = lse(prob_cond) - log(B)
      # average = exp(average_log)
      # marglik = marglik * average
      marglik_log = marglik_log + average_log
    }
    # weights = exp(prob_cond)/(sum(exp(prob_cond)))
    weights_log = prob_cond - lse(prob_cond)
    weights = exp(weights_log)
    
    if(all(is.na(weights))){
      # zero weights for every particles, then sample with equal weights
      label = sample(1:B, B, replace = T)
      cat('Caution: all particle weights are zero!\n')
    }else{
      label = sample(1:B, B, replace = T, prob = weights)
    }
    
    # generate particles at m
    if(parallel){
      clusterSize = as.numeric(Sys.getenv("NSLOTS"))
      if (is.na(clusterSize)) clusterSize = 1
      cluster = makeCluster(clusterSize, type = "SOCK")
      clusterExport(cluster, "N")
      results = clusterApply(cluster, 1:clusterSize, nextBpart_splited, 
                             particles=particles, label=label, W=W, t_obs=t_obs, p=para_est[1], q=para_est[2], 
                             gamma=para_est[3], B=B, m=m, N=N, edgeIndex=edgeIndex, process, clusterSize, seed_parallel)
      portion = B/clusterSize
      for (cl in 1:clusterSize){
        for(b in 1:portion){
          particles[[(cl-1)*portion+b]][[m]] = results[[cl]]$part[[b]]
          W[(cl-1)*portion+b, m] = results[[cl]]$W[b]
        }
      }
      stopCluster(cluster)
    }else{
      for (b in 1:B){
        #net_start based on weights
        net_start = particles[[label[b]]][[m-1]]
        w_start = W[label[b], m-1]
        temp = get_SimSeq2(para_est[1], para_est[2], para_est[3], N, t_obs[m]-t_obs[m-1], net_start, w_start, process = process, edgeIndex = edgeIndex)
        particles[[b]][[m]] = temp$network_seq_obs[[1]]
        # make sure None is appended 
        #if (length(temp$network_seq_obs[[1]]) == 0) particles[[b]][[m]] = integer(0)
        W[b, m] = temp$bin_var_obs
      }
    }
    
    # delete particles at generation m-1
    for (index in 1:B){
      particles[[index]][[m-1]] = '*'
    }
  }
  
  # average at time m
  prob_cond = sapply(particles, FUN = function(list){
    res = getCondProb(observe = Y[[M]], true = format_back(list[[M]], N), para_est[4], para_est[5])
    return (res)})
  average_log = lse(prob_cond) - log(B)
  # average = exp(average_log)
  # marglik = marglik * average
  marglik_log = marglik_log + average_log
  return (marglik_log)
}
  