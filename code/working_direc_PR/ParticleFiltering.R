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
particleFil = function (B, Y, t_obs, p, q, gamma, alpha, beta, process = "ER", E_step = T){
  N = nrow(Y[[1]])
  M = length(t_obs)
  W = matrix(0, B, M) # each row is the seq of binary var at t_obs
  W[,1] = rep(1, B) # set binary var to be 1 at the first obs. moments
  
  # create B particles at each t_obs, store sequential particles at t_obs in a list 
  # particles is a list of list, particles[[1]],...,particles[[B]]
  particles = array(list(), B)
  for (b in 1:B){
    particles[[b]] = list()
    particles[[b]][[1]] = Y[[1]]
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
      res = getCondProb(observe = Y[[m-1]], true = list[[m-1]], alpha, beta) 
      return (res)})
    if (E_step){
      weights = prob_cond/(sum(prob_cond))
    }
    else{
      weights = exp(prob_cond)/(sum(exp(prob_cond)))
      # weights_log = prob_cond - lse(prob_cond)
      # weights = exp(weights_log)
    }
    #sorted = sort(weights, decreasing = T, index.return = T)
    #label = rep(sorted$ix, round(sorted$x*B,digits = 0))
    label = sample(1:B, B, replace = T, prob = weights)
    ancesters[, m] = label
    for (b in 1:B){
      #net_start based on weights
      net_start = particles[[label[b]]][[m-1]]
      w_start = W[label[b], m-1]
      temp = get_SimSeq(p, q, gamma, N, t_obs[m]-t_obs[m-1], net_start, w_start, process = process)
      particles[[b]][[m]] = temp$network_seq_obs[[1]]
      W[b, m] = temp$bin_var_obs
      # num_tran = length(temp$t_transition)
      # W[b,m] = checkAdding(temp$network_seq_hidden[[num_tran-1]],
      #                      temp$network_seq_hidden[[num_tran]])
    }
  }
  prob_cond = sapply(particles, FUN = function(list){
    res = getCondProb(observe = Y[[M]], true = list[[M]], alpha, beta)
    return (res)})
  weights = prob_cond/(sum(prob_cond))
  return (list(particles = particles, W = W, weights = weights, ancesters = ancesters))
}

# particleFil = function (B, Y, t_obs, p, q, gamma, alpha, beta, process = "ER"){
#   N = nrow(Y[[1]])
#   M = length(t_obs)
#   W = matrix(0, B, M) # each row is the seq of binary var at t_obs
#   W[,1] = rep(1, B) # set binary var to be 1 at the first obs. moments
#   
#   # create B particles at each t_obs, store sequential particles at t_obs in a list 
#   # particles is a list of list, particles[[1]],...,particles[[B]]
#   particles = array(list(), B)
#   for (b in 1:B){
#     particles[[b]] = list()
#     particles[[b]][[1]] = Y[[1]]
#   }
#   ancesters = matrix(0, B, M)
#   ancesters[, 1] = 1:B
#   for (m in 2:M){
#     cat("m = ", m, "...")
#     # if (m%%20 == 0){
#     #   cat("m = ", m, "...")
#     # }
#     # check weights
#     prob_cond = sapply(particles, FUN = function(list){
#       res = getCondProb(observe = Y[[m-1]], true = list[[m-1]], alpha, beta) 
#       return (res)})
#     weights = prob_cond/(sum(prob_cond))
#     # weights = exp(prob_cond)/(sum(exp(prob_cond)))
#     # weights_log = prob_cond - lse(prob_cond)
#     # weights = exp(weights_log)
#     #sorted = sort(weights, decreasing = T, index.return = T)
#     #label = rep(sorted$ix, round(sorted$x*B,digits = 0))
#     label = sample(1:B, B, replace = T, prob = weights)
#     ancesters[, m] = label
#     for (b in 1:B){
#       #net_start based on weights
#       net_start = particles[[label[b]]][[m-1]]
#       w_start = W[label[b], m-1]
#       temp = get_SimSeq(p, q, gamma, N, t_obs[m]-t_obs[m-1], net_start, w_start, process = process)
#       particles[[b]][[m]] = temp$network_seq_obs[[1]]
#       W[b, m] = temp$bin_var_obs
#       # num_tran = length(temp$t_transition)
#       # W[b,m] = checkAdding(temp$network_seq_hidden[[num_tran-1]],
#       #                      temp$network_seq_hidden[[num_tran]])
#     }
#   }
#   prob_cond = sapply(particles, FUN = function(list){
#     res = getCondProb(observe = Y[[M]], true = list[[M]], alpha, beta)
#     return (res)})
#   weights = prob_cond/(sum(prob_cond))
#   return (list(particles = particles, W = W, weights = weights, ancesters = ancesters))
# }

# ancestral sampling
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
# t_obs
# Y = net_obs_noise
# particlesRes = particleFil (B = 3, Y, t_obs, p, q, gamma, alpha, beta,
#                        process = "ER")
# particles = particlesRes$particles
# W=particlesRes$W
# particles[[1]]
# W

#wrapper funciton of particleFil for parallelization
particleFil_parallel = function (B, Y, t_obs, p, q, gamma, alpha, beta, process = "ER"){
  #theta = c(p, q, gamma, alpha, beta)
  N = nrow(Y[[1]])
  M = length(t_obs)
  W = matrix(0, B, M) #each row is the seq of binary var at t_obs
  W[,1] = rep(1, B) #set binary var to be 1 at the first obs. moments
  
  #create B particles at each t_obs, store sequential particles at t_obs in a list 
  #particles is a list of list, particles[[1]],...,particles[[B]]
  particles = array(list(), B)
  for (b in 1:B){
    particles[[b]] = list()
    particles[[b]][[1]] = Y[[1]]
  }
  # cat("m=2...")
  for (m in 2:M){
    cat("m = ", m, "...")
    # if (m%%20 == 0){
    #   cat("m = ", m, "...")
    # }
    #check weights
    prob_cond = sapply(particles, FUN = function(list){
      res = getCondProb (observe = Y[[m-1]], true = list[[m-1]], alpha, beta) 
      return (res)})
    weights = prob_cond/(sum(prob_cond))
    #sorted = sort(weights, decreasing = T, index.return = T)
    #label = rep(sorted$ix, round(sorted$x*B,digits = 0))
    label = sample(1:B, B, replace = T, prob = weights)
    clusterSize = as.numeric(Sys.getenv("NSLOTS"))
    if (is.na(clusterSize)) clusterSize = 1
    cluster = makeCluster(clusterSize, type = "SOCK")
    results = clusterApply(cluster, 1:B, nextBpart_splited, 
                           particles, label, W, t_obs, p, q, gamma, B, m, N)
    #portion = B/clusterSize
    for (cl in 1:B){
      particles[[cl]][[m]] = results[[cl]]$part
      W[cl, m] = results[[cl]]$W
    }
    stopCluster(cluster)
  }
  return (list(particles=particles, W = W))
}

nextBpart_splited = function(processID, particles, label, W, t_obs, p, q, gamma, B, m, N){
  source("GetSimSeq.R")
  part = list()
  WW= 0 
  id = 1
  b = processID
  #net_start based on weights
  net_start=particles[[label[b]]][[m-1]]
  w_start=W[label[b],m-1]
  temp = get_SimSeq(p, q, gamma, N, t_obs[m]-t_obs[m-1], net_start, w_start)
  part = temp$network_seq_obs[[1]]
  WW = temp$bin_var_obs
  # num_tran = length(temp$t_transition)
  # WW = checkAdding(temp$network_seq_hidden[[num_tran-1]],       
  #                   temp$network_seq_hidden[[num_tran]])
  return (list(part = part, W=WW))
} 

