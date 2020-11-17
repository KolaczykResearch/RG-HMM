source("utilities.R")
source("prob_fun.R")
source("GetSimSeq.R")
  
# particleFil()----------------------------------
# Get B particles at each observational moments by ER
#-------Input----------------------------------
# B: number of particles at each observational moments
# Y: observed sequence of networks (with noise)
# t_obs: a vetor of observational time points
# p: birth rate
# q: death rate 
# gamma: rate parameter 
# alpha: birth rate 
# beta: death rate
# process: "ER" or "PR"
# (M: total number of observational time points)
# (N: total number of vertices)
#------Output---------------------------------
# particles: an array of B list, each list stores a sequence of true networks at observational moments
# W: B by M matrix, each row represents a sequence of binary variables at M observational moments
particleFil = function (B, Y, t_obs, p, q, gamma, alpha, beta, process = "ER", edgeIndex){
  N = nrow(Y[[1]])
  M = length(t_obs)
  W = matrix(0, B, M) # each row is the seq of binary var at t_obs
  if (is.empty(Y[[1]])){
    W[,1] = rep(0, B)
  }else{
    W[,1] = rep(1, B) # set binary var to be 1 at the first obs. moments if Y[[1]] is not empty, else 0
  }
  
  # create B particles at each t_obs, store sequential particles at t_obs in a list 
  # particles is a list of list, particles[[1]],...,particles[[B]]
  particles = array(list(), B)
  Y1 = reformat(Y[[1]], N)
  for (b in 1:B){
    particles[[b]] = list()
    # particles[[b]][[1]] = Y[[1]]
    particles[[b]][[1]] = Y1
  }
  ancesters = matrix(0, B, M)
  ancesters[, 1] = 1:B
  for (m in 2:M){
    cat("m = ", m, "...")
    # if (m%%20 == 0){
    #   cat("m = ", m, "...")
    # }
    # check weights
    prob_cond = sapply(particles, FUN = function(list){
      res = getCondProb(observe = Y[[m-1]], true = format_back(list[[m-1]], N), alpha, beta) 
      return (res)})
    # weights = exp(prob_cond)/(sum(exp(prob_cond)))
    weights_log = prob_cond - lse(prob_cond)
    weights = exp(weights_log)
    #sorted = sort(weights, decreasing = T, index.return = T)
    #label = rep(sorted$ix, round(sorted$x*B,digits = 0))
    #cat(weights, '\n')
    label = sample(1:B, B, replace = T, prob = weights)
    ancesters[, m] = label
    
    # delete particles with no offsprings, as well as their ancesters
    rmset = setdiff(1:B, unique(label))
    for (index in rmset){
      particles[[index]][[m-1]] = '*'
    }
    m_ = m
    while (m_ > 2){
      rmset = setdiff(ancesters[rmset, m_-1], ancesters[-rmset, m_-1])
      for (index in rmset){
        particles[[index]][[m_-2]] = '*'
      }
      m_ = m_ - 1
    }

    for (b in 1:B){
      #net_start based on weights
      net_start = particles[[label[b]]][[m-1]]
      w_start = W[label[b], m-1]
      #temp = get_SimSeq(p, q, gamma, N, t_obs[m]-t_obs[m-1], format_back(net_start, N), w_start, process = process)
      temp = get_SimSeq2(p, q, gamma, N, t_obs[m]-t_obs[m-1], net_start, w_start, process = process, edgeIndex = edgeIndex)
      particles[[b]][[m]] = temp$network_seq_obs[[1]]
      W[b, m] = temp$bin_var_obs
      # num_tran = length(temp$t_transition)
      # W[b,m] = checkAdding(temp$network_seq_hidden[[num_tran-1]],
      #                      temp$network_seq_hidden[[num_tran]])
    }
  }
  prob_cond = sapply(particles, FUN = function(list){
    res = getCondProb(observe = Y[[M]], true = format_back(list[[M]], N), alpha, beta)
    return (res)})
  # weights = prob_cond/(sum(prob_cond))
  # weights = exp(prob_cond)/(sum(exp(prob_cond)))
  weights_log = prob_cond - lse(prob_cond)
  weights = exp(weights_log)
  return (list(particles = particles, W = W, weights = weights, ancesters = ancesters))
}

# ancestral sampling (Backward Sampling)
# given the index of the last particle, track the ancesters
# return a particle path
get_particle_path = function(partic, index, M){
  Y_true = array(list(), M)
  W_true = numeric(M)
  Y_true[[M]] = partic$particles[[index]][[M]]
  W_true[M] = partic$W[index, M]
  for (m in (M-1):1){
    index = partic$ancesters[index, m+1]
    Y_true[[m]] = partic$particles[[index]][[m]]
    W_true[m] = partic$W[index, m]
  }
  return (list(Y_true = Y_true, W_true=W_true))
}

# Wrapper funciton of particleFil with parallelization
particleFil_parallel = function (B, Y, t_obs, p, q, gamma, alpha, beta, process = "ER", edgeIndex, seed_parallel=NULL){
  N = nrow(Y[[1]])
  M = length(t_obs)
  W = matrix(0, B, M) # each row is the seq of binary var at t_obs
  if (is.empty(Y[[1]])){
    W[,1] = rep(0, B)
  }else{
    W[,1] = rep(1, B) # set binary var to be 1 at the first obs. moments if Y[[1]] is not empty, else 0
  }
  
  # create B particles at each t_obs, store sequential particles at t_obs in a list 
  # particles is a list of list, particles[[1]],...,particles[[B]]
  particles = array(list(), B)
  Y1 = reformat(Y[[1]], N)
  for (b in 1:B){
    particles[[b]] = list()
    particles[[b]][[1]] = Y1
  }
  ancesters = matrix(0, B, M)
  ancesters[, 1] = 1:B
  for (m in 2:M){
    cat("m = ", m, "...")

    # check weights and resampling
    prob_cond = sapply(particles, FUN = function(list){
      res = getCondProb(observe = Y[[m-1]], true = format_back(list[[m-1]], N), alpha, beta) 
      return (res)})
    weights_log = prob_cond - lse(prob_cond)
    weights = exp(weights_log)
    label = sample(1:B, B, replace = T, prob = weights)
    ancesters[, m] = label
    
    # delete particles with no offsprings, as well as their ancesters
    rmset = setdiff(1:B, unique(label))
    for (index in rmset){
      particles[[index]][[m-1]] = '*'
    }
    m_ = m
    while (m_ > 2){
      rmset = setdiff(ancesters[rmset, m_-1], ancesters[-rmset, m_-1])
      for (index in rmset){
        particles[[index]][[m_-2]] = '*'
      }
      m_ = m_ - 1
    }
    
    # propagate
    clusterSize = as.numeric(Sys.getenv("NSLOTS"))
    if (is.na(clusterSize)) clusterSize = 1
    cluster = makeCluster(clusterSize, type = "SOCK")
    clusterExport(cluster, "N")
    results = clusterApply(cluster, 1:clusterSize, nextBpart_splited, 
                           particles=particles, label=label, W=W, t_obs=t_obs, p=p, q=q, 
                           gamma=gamma, B=B, m=m, N=N, edgeIndex=edgeIndex, process, clusterSize, seed_parallel)
    portion = B/clusterSize
    for (cl in 1:clusterSize){
      for(b in 1:portion){
        particles[[(cl-1)*portion+b]][[m]] = results[[cl]]$part[[b]]
        W[(cl-1)*portion+b, m] = results[[cl]]$W[b]
      }
    }
    stopCluster(cluster)
  }
  prob_cond = sapply(particles, FUN = function(list){
    res = getCondProb(observe = Y[[M]], true = format_back(list[[M]], N), alpha, beta)
    return (res)})
  weights_log = prob_cond - lse(prob_cond)
  weights = exp(weights_log)
  return (list(particles = particles, W = W, weights = weights, ancesters = ancesters))
}

nextBpart_splited = function(processID, particles, label, W, t_obs, p, q, gamma, B, m, N, 
                             edgeIndex, process, clusterSize, seed_parallel=NULL){
  source("Model/GetSimSeq.R", chdir=T)
  # set seed for the env if seed is given
  if(length(seed_parallel)!=0){
    set.seed(seed_parallel+processID)
  }
  portion = B/clusterSize
  part = list()
  WW= 0 
  id = 1
  for (b in ((processID-1)*portion+1):(processID*portion)){
    #net_start based on weights
    net_start = particles[[label[b]]][[m-1]]
    w_start = W[label[b], m-1]
    temp = get_SimSeq2(p, q, gamma, N, t_obs[m]-t_obs[m-1], net_start, w_start, process = process, edgeIndex = edgeIndex)
    part[[id]] = temp$network_seq_obs[[1]]
    WW[id] = temp$bin_var_obs
    id = id+1
  }
  return (list(part = part, W=WW))
} 

