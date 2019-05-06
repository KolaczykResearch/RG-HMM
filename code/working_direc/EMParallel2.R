source("ParticleFiltering.R")
source("MHSamples2.R")
library("snow")

#particles: B sampled true network in the hidden layer
#Y: observed network
#burnIn: the number of first few obs. moments at which particles are excluded from updating alpha, beta
#phi: the number of true network sample series used for updating alpha, beta
# update_alpha_beta = function(particles, Y, burnIn, phi){
#   B = length(particles)
#   M = length(particles[[1]])
#   N = nrow(particles[[1]][[1]])
#   set = sample(B, phi)
#   sum_c = 0
#   sum_cd = 0
#   sum_b = 0
#   sum_ab = 0
#   for (i in set){
#     for (j in (burnIn+1):M){
#       sum_c = sum_c + sum((Y[[j]] - particles[[i]][[j]])==1)/2
#       sum_cd = sum_cd + nonedge_num(particles[[i]][[j]])
#       sum_b = sum_b + sum((Y[[j]] - particles[[i]][[j]])==-1)/2
#       sum_ab = sum_ab + edge_num(particles[[i]][[j]])
#     }
#   }
#   alpha_hat = sum_c/sum_cd
#   beta_hat = sum_b/sum_ab
#   cat("alpha = ", alpha_hat, ", beta = ", beta_hat)
#   return (c(alpha_hat, beta_hat))
# }

update_alpha_beta = function(partic, Y, burnIn, phi){
  B = length(partic$particles)
  M = length(partic$particles[[1]])
  N = nrow(partic$particles[[1]][[1]])
  # set = sample(B, phi)
  # sample phi last indices
  set = sample(1:B, phi, replace = F, prob = partic$weights) 
  sum_c = 0
  sum_cd = 0
  sum_b = 0
  sum_ab = 0
  for (i in set){
    temp = get_particle_path(partic, i, M)
    Y_true = temp$Y_true
    for (j in (burnIn+1):M){
      sum_c = sum_c + sum((Y[[j]] - Y_true[[j]])==1)/2
      sum_cd = sum_cd + nonedge_num(Y_true[[j]])
      sum_b = sum_b + sum((Y[[j]] - Y_true[[j]])==-1)/2
      sum_ab = sum_ab + edge_num(Y_true[[j]])
    }
  }
  alpha_hat = sum_c/sum_cd
  beta_hat = sum_b/sum_ab
  cat("alpha = ", alpha_hat, ", beta = ", beta_hat)
  return (c(alpha_hat, beta_hat))
}

update_pqgamma = function(counts_vec, H, D, t_obs, M){
  T = t_obs[M] - t_obs[1]
  sum_01 = counts_vec[1]
  sum_0 = counts_vec[2]
  sum_10= counts_vec[3]
  sum_1 = counts_vec[4] 
  sum_R = counts_vec[5]
  p_hat = sum_01 / sum_0
  q_hat = sum_10 / sum_1
  gamma_hat = (sum_R/(H*D))/T
  cat("p = ", p_hat, ", q = ", q_hat, ", gamma = ", gamma_hat, "\n")
  return (c(p_hat, q_hat, gamma_hat))
}

# EM()-----------------------------------
# -------Input----------------------------------
# Y: Observed sequence of networks (with noise)
# t_obs: A vetor of observational time points
# B: Number of particles generated at each observational moments
# (M: Total number of observational time points)
# (N: Total number of vertices)
# process: "ER" or "PR"
# thr: Stopping threshold for EM iteration
# init: Initialized value for p, q, gamma, alpha, beta
# D: The length of the chain for mcmc 
# MaxIter: Maximum number of iterations allowed for EM 
# prop_max: accept proposed sample path when prop_max exceeded
# burnIn: Number of observational moments dropped to compute alpha and gamma
# burnIn_mcmc: Number of observational moments dropped in mcmc
# phi: Number of particles used to compute alpha and gamma
# H: Number of true network sequences used in mcmc 
# particleParallel: TRUE if running particle filtering in parallel
# MHParallel: TRUE if running mcmc in parallel

#step_init: 
#n1: number of sample paths used to estimate the derivative of expectation
#n2: the number of updates of parameters before averaging
#times:
#MaxIter2: 
#clusterSizePF: number of cores for particle filtering (exact division of M)
#clusterSizeRM: number of cores for RM (exact division of M)
#outFile: name of outFile
# ------Output---------------------------------
# estimation: Vector of p, q, gamma, alpha, beta
# iteration: Number of EM iterations used for convergence
EM_parallel = function(Y, t_obs, B, process = "ER", thr=1e-1, init = NULL, D = 50, MaxIter = 300,
                       prop_max = 100, burnIn = 50, phi, H = 100,
                       particleParallel = FALSE, MHParallel = FALSE, burnIn_mcmc = 0){
  if (length(init)==0){
    init = rep(0.1,5)
  }
  p_cur = init[1]
  q_cur = init[2]
  gamma_cur = init[3]
  alpha_cur = init[4]
  beta_cur = init[5]
  cat("Initialization:", init, "\n")
  M = length(t_obs)
  # stopping condition for EM:
  stopConEM = 1
  iter = 0
  
  #keep track of each EM updates??? ####################################
  while (stopConEM > thr & iter < MaxIter){
    # Particle Filtering
    cat("---------------------------\n")
    cat("EM Iteration",iter+1,": \n--> Generating particles...\n")
    time1 = proc.time()
    if(particleParallel){
      partic = particleFil_parallel (B, Y, t_obs, p_cur, q_cur, gamma_cur, alpha_cur, beta_cur,
                                     process = process)
    }else{
      partic = particleFil(B, Y, t_obs, p_cur, q_cur, gamma_cur, alpha_cur, beta_cur,
                           process = process)
    }
    cat("\n  Finished particle filtering. \n")
    time_used2 = proc.time() - time1
    cat("  Time elapsed:", time_used2[3],"\n")
    
    # E-step (MCMC)
    cat("--> Start E-step...\n")
    time1 = proc.time()
    # Sample H true network in the hidden layer from B particles
    # then for each of them sample D paths that connects the truth
    # \hat p = E[E[\sum(I{w_{n-1} = 0, w_m = 1)]]/E[E[\sum(I_{w_{n-1} = 0})]]
    list_sum = list()
    #pqgamma_track = matrix(0, H, 3)
    counts_vec = rep(0, 5)
    last_particle_index = sample(1:B, H, replace = F, prob = partic$weights)
    cat("   Sampling MH paths for b = ...\n   ")
    for (h in 1:H){
      index = last_particle_index[h]
      # for checking purose
      list_sum[[h]] = matrix(0, M-1-burnIn, 5)
      cat(h, ", ")
      # cat(h, ":", ids[h], "...,")
      # b = ids[h]
      # Y_true = partic$particles[[b]]
      # W_true = partic$W[b,]
      temp = get_particle_path(partic, index, M)
      Y_true = temp$Y_true
      W_true = temp$W_true
      # sample D paths that connects the truth, (parallel for M here, not D)
      if (MHParallel){
        clusterSize = as.numeric(Sys.getenv("NSLOTS"))
        if (is.na(clusterSize)) clusterSize = 1
        #initialize a cluster for parallel computing
        cluster = makeCluster(clusterSize, type = "SOCK")
        #sample D paths that connects Y_true, W_true
        results = clusterApply(cluster, 1:(M-1-burnIn), getSamplePath_parallel, D, Y_true, W_true, t_obs,
                               p_cur, q_cur, gamma_cur, process = process, prop_max, burnIn, burnIn_mcmc)
        # counts_vec = c(sum_01, sum_0, sum_10, sum_1, sum_R)
        for (cl in 1:(M-1-burnIn)){
          counts_vec = counts_vec + results[[cl]]
          list_sum[[h]][cl, ] = results[[cl]]
        }
        stopCluster(cluster)
      }
      else{
        #sample D paths that connects Y_true, W_true
        for (m in (burnIn+1):(M-1)){
          results = getSamplePath_nonparallel(m, D, Y_true, W_true, t_obs,
                                           p_cur, q_cur, gamma_cur, process = process, prop_max, burnIn_mcmc)
          counts_vec = counts_vec + results
          list_sum[[h]][m-burnIn, ] = results
        }
      }
     # pqgamma_track[h, ] = update_pqgamma(counts_vec, H=h, D, t_obs, M)
    }
    time_used2 = proc.time() - time1
    cat("  Finished E-step.\n")
    cat("  Time elapsed:", time_used2[3],"\n")
    
    # M-Step
    #EM update
    cat("--> Start M-step...\n")

    cat('\n Updating p, q, gamma... \n')
    pqgamma = update_pqgamma(counts_vec, H, D, t_obs, M)
    
    cat('\n Updating alpha, beta... \n')
    ##
    #partic = particleFil(B, Y, t_obs, p_cur, q_cur, gamma_cur, alpha_cur, beta_cur,
    #                     process = process, E_step = F)
    time1 = proc.time()
    alphaBeta = update_alpha_beta(partic, Y, burnIn, phi)
    cat("\n  Finished updating. \n")
    time_used2 = proc.time() - time1
    cat("  Time elapsed:", time_used2[3],"\n")
  
    p_nex = pqgamma[1]
    q_nex = pqgamma[2]
    gamma_nex = pqgamma[3]
    alpha_nex = alphaBeta[1]
    beta_nex = alphaBeta[2]
    
    # Check EM stopping condition
    stopConEM = sqrt(sum((c(p_nex, q_nex, gamma_nex, alpha_nex,
                            beta_nex)-c(p_cur, q_cur, gamma_cur, alpha_cur, beta_cur))^2))/sqrt(sum((c(p_cur, q_cur, gamma_cur, alpha_cur, beta_cur))^2))
    cat("\n Relative change:", stopConEM, "\n")
    p_cur = p_nex
    q_cur = q_nex
    gamma_cur = gamma_nex
    alpha_cur = alpha_nex
    beta_cur = beta_nex
    iter = iter + 1
  }
  cat("---------------------------\n")
  cat("EM Iterations Done!\n")
  return (list(estimation = c(p_cur, q_cur, gamma_cur, alpha_cur, beta_cur),
               iteration = iter, list_sum = list_sum))
}


# EM_parallel = function(Y, t_obs, B, process = "ER", thr=1e-1, init = NULL, D = 50, MaxIter = 300, 
#                        prop_max = 100, burnIn = 50, phi, H = 100, 
#                        particleParallel = FALSE, MHParallel = FALSE, burnIn_mcmc = 0){
#   if (length(init)==0){
#     init = rep(0.1,5)
#   }
#   p_cur = init[1]
#   q_cur = init[2]
#   gamma_cur = init[3]
#   alpha_cur = init[4]
#   beta_cur = init[5]
#   cat("Initialization:", init, "\n")
#   M = length(t_obs)
#   # stopping condition for EM:
#   stopConEM = 1
#   iter = 0
#   
#   #keep track of each EM updates??? ####################################
#   while (stopConEM > thr & iter < MaxIter){
#     # Particle Filtering
#     cat("---------------------------\n")
#     cat("EM Iteration",iter+1,": \n--> Generating particles...\n")
#     time1 = proc.time()
#     if(particleParallel){
#       partic = particleFil_parallel (B, Y, t_obs, p_cur, q_cur, gamma_cur, alpha_cur, beta_cur, 
#                                      process = process)
#     }else{
#       partic = particleFil(B, Y, t_obs, p_cur, q_cur, gamma_cur, alpha_cur, beta_cur, 
#                            process = process)
#     }
#     cat("\n  Finished particle filtering. \n")
#     time_used2 = proc.time() - time1
#     cat("  Time elapsed:", time_used2[3],"\n")
#     
#     # E-step (MCMC)
#     cat("--> Start E-step...\n")
#     time1 = proc.time()
#     # Sample H true network in the hidden layer from B particles
#     # then for each of them sample D paths that connects the truth
#     # \hat p = E[E[\sum(I{w_{n-1} = 0, w_m = 1)]]/E[E[\sum(I_{w_{n-1} = 0})]]
#     counts_vec = rep(0, 5)
#     ids = sample(1:B, H)
#     if(MHParallel){
#       cat("   Sampling MH paths for b = ...\n   ")
#       clusterSize = as.numeric(Sys.getenv("NSLOTS"))
#       if (is.na(clusterSize)) clusterSize = 1
#       #initialize a cluster for parallel computing
#       cluster = makeCluster(clusterSize, type = "SOCK")
#       #sample D paths that connects Y_true, W_true
#       results = clusterApply(cluster, 1:(H*(M-1-burnIn)), getSamplePath_parallel, D, partic, t_obs,
#                              p_cur, q_cur, gamma_cur, process = process, prop_max, ids, burnIn, burnIn_mcmc)
#       # counts_vec = c(sum_01, sum_0, sum_10, sum_1, sum_R)
#       for (cl in 1:(H*(M-1-burnIn))){
#         counts_vec = counts_vec + results[[cl]]
#       }
#       stopCluster(cluster)
#     }else{
#       cat("   Sampling MH paths for b = ...\n   ")
#       for (h in 1:H){
#         cat(h, ":", ids[h], "...,")
#         b = ids[h]
#         Y_true = partic$particles[[b]]
#         W_true = partic$W[b, ]
#         #sample D paths that connects Y_true, W_true
#         for (m in (burnIn+1):(M-1)){
#           results = getSamplePath_nonparallel(m, D, Y_true, W_true, t_obs,
#                                            p_cur, q_cur, gamma_cur, process = process, prop_max, burnIn_mcmc)
#           counts_vec = counts_vec + results
#         }
#       }
#     }
#     time_used2 = proc.time() - time1
#     cat("  Finished E-step.\n")
#     cat("  Time elapsed:", time_used2[3],"\n")
#     
#     # M-Step
#     #EM update
#     cat("--> Start M-step...\n")
#     pqgamma = update_pqgamma(counts_vec, H, D, t_obs, M)
#     alphaBeta = update_alpha_beta(partic$particles, Y, burnIn, phi)
#     p_nex = pqgamma[1]
#     q_nex = pqgamma[2]
#     gamma_nex = pqgamma[3]
#     alpha_nex = alphaBeta[1]
#     beta_nex = alphaBeta[2]
#     
#     # Check EM stopping condition
#     stopConEM = sqrt(sum((c(p_nex, q_nex, gamma_nex, alpha_nex,
#                             beta_nex)-c(p_cur, q_cur, gamma_cur, alpha_cur, beta_cur))^2))/sqrt(sum((c(p_cur, q_cur, gamma_cur, alpha_cur, beta_cur))^2))
#     cat("\n Relative change:", stopConEM, "\n")
#     p_cur = p_nex
#     q_cur = q_nex
#     gamma_cur = gamma_nex
#     alpha_cur = alpha_nex
#     beta_cur = beta_nex
#     iter = iter + 1
#   }
#   cat("---------------------------\n")
#   cat("EM Iterations Done!\n")
#   return (list(estimation = c(p_cur, q_cur, gamma_cur, alpha_cur, beta_cur),
#                iteration = iter))
# }





