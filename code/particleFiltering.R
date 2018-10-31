source("utilities.R")
source("probFunc.R")
source("GetSimSeq.R")

#particleFil()-----------------------------------
#-------Input----------------------------------
#B: number of particles at each observational moments
#Y: observed sequence of networks (with noise)
#t_obs: a vetor of observational time points
#p: birth rate
#q: death rate 
#gamma: rate parameter 
#alpha: birth rate 
#beta: death rate
#process: "ER" or "PR"
#(M: total number of observational time points)
#(N: total number of vertices)
#------Output---------------------------------
##get B particles at each observation moments by ER
#Y: observed seq of networks
particleFil = function (B, Y, t_obs, p, q, gamma, alpha, beta, 
                        process = "ER", outFile = "output.txt"){
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
  cat("m=2...")
  cat("m=2...", sep = " ", file = outFile, append = TRUE)
  for (m in 2:M){
    cat("m = ", m, "...")
    # if (m%%20 == 0){
    #   cat("m = ", m, "...")
    #   cat("m = ", m, "...", sep = " ", file = outFile, append = TRUE)
    # }
    #check weights
    prob_cond = sapply(particles, FUN = function(list){
      res = getCondProb (observe = Y[[m-1]], true = list[[m-1]], alpha, beta) 
      return (res)})
    weights = prob_cond/(sum(prob_cond))
    sorted = sort(weights, decreasing = T, index.return = T)
    label = rep(sorted$ix, round(sorted$x*B,digits = 0))
    for (b in 1:B){
      #net_start based on weights
      #cat(label[b], " ")
      net_start=particles[[label[b]]][[m-1]]
      w_start=W[label[b],m-1]
      temp = get_SimSeq(p, q, gamma, N, t_obs[m]-t_obs[m-1], net_start, w_start)
      particles[[b]][[m]] = temp$network_seq_obs[[1]]
      num_tran = length(temp$t_transition)
      W[b,m] = checkAdding(temp$network_seq_hidden[[num_tran-1]], 
                           temp$network_seq_hidden[[num_tran]])
    }
    #cat("\n")
  }
  return (list(particles=particles, W = W))
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
particleFil_parallel = function (B, Y, t_obs, p, q, gamma, alpha, beta, 
                                 process = "ER", outFile = "output.txt"){
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
  cat("m=2...")
  for (m in 2:M){
    if (m%%20 == 0){
      cat("m = ", m, "...")
      cat("m = ", m, "...", sep = "\n", file = outFile, append = TRUE)
    }
    #check weights
    prob_cond = sapply(particles, FUN = function(list){
      res = getCondProb (observe = Y[[m-1]], true = list[[m-1]], alpha, beta) 
      return (res)})
    weights = prob_cond/(sum(prob_cond))
    sorted = sort(weights, decreasing = T, index.return = T)
    label = rep(sorted$ix, round(sorted$x*B,digits = 0))
    
    clusterSize = as.numeric(Sys.getenv("NSLOTS"))
    if (is.na(clusterSize)) clusterSize = 1

    cluster = makeCluster(clusterSize, type = "SOCK")
    # results = clusterApply(cluster, 1:B, nextBpart_splited, 
    #                        particles, label, W, t_obs, p, q, gamma, B, m, N, clusterSize)
    results = clusterApply(cluster, 1:B, nextBpart_splited, 
                           particles, label, W, t_obs, p, q, gamma, B, m, N)
    #portion = B/clusterSize
    for (cl in 1:B){
        particles[[cl]][[m]] = results[[cl]]$part
        W[cl, m] = results[[cl]]$W
      #for(b in 1:portion){
       # particles[[(cl-1)*portion+b]][[m]] = results[[cl]]$part[[b]]
       # W[(cl-1)*portion+b, m] = results[[cl]]$W[b]
      #}
    }
    stopCluster(cluster)
  }
  return (list(particles=particles, W = W))
}

nextBpart_splited = function(processID, particles, label, W, t_obs, p, q, gamma, B, m, N){
  source("GetSimSeq.R")
  #portion = B/clusterSize
  part = list()
  WW= 0 
  id = 1
  b = processID
  #net_start based on weights
  #cat(label[b], " ")
  net_start=particles[[label[b]]][[m-1]]
  w_start=W[label[b],m-1]
  temp = get_SimSeq(p, q, gamma, N, t_obs[m]-t_obs[m-1], net_start, w_start)
  part = temp$network_seq_obs[[1]]
  num_tran = length(temp$t_transition)
  WW = checkAdding(temp$network_seq_hidden[[num_tran-1]],       
                    temp$network_seq_hidden[[num_tran]])
  return (list(part = part, W=WW))
} 

