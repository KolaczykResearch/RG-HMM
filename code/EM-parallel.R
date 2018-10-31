source("particleFiltering.R")
source("MHSamples.R")
# source("RobbinsMonro.R")
library("snow")

#particles: B sampled true network in the hidden layer
#Y: observed network
#burnIn: the number of first few obs. moments at which particles are excluded from updating alpha, beta
#phi: the number of true network sample series used for updating alpha, beta
update_alpha_beta = function(particles, Y, burnIn, phi){
  B = length(particles)
  M = length(particles[[1]])
  N = nrow(particles[[1]][[1]])
  set = sample(B, phi)
  sum_c = 0
  sum_cd = 0
  sum_b = 0
  sum_ab = 0
  for (i in set){
    for (j in burnIn:M){
      sum_c = sum_c + sum((Y[[j]] - particles[[i]][[j]])==1)/2
      sum_cd = sum_cd + nonedge_num(particles[[i]][[j]])
      sum_b = sum_b + sum((Y[[j]] - particles[[i]][[j]])==-1)/2
      sum_ab = sum_ab + edge_num(particles[[i]][[j]])
    }
  }
  alpha_hat = sum_c/sum_cd
  beta_hat = sum_b/sum_ab
  cat("alpha = ", alpha_hat, ", beta = ", beta_hat)
  cat("alpha = ", alpha_hat, ", beta = ", beta_hat, sep = " ", file = outFile, append = TRUE)
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
  cat("p = ", p_hat, ", q = ", q_hat, "gamma = ", gamma_hat, "\n")
  cat("p = ", p_hat, ", q = ", q_hat, "gamma = ", gamma_hat, "\n", file = outFile, append = TRUE)
  return (c(p_hat, q_hat, gamma_hat))
}

#EM()-----------------------------------
#-------Input----------------------------------
#Y: observed sequence of networks (with noise)
#t_obs: a vetor of observational time points
#B: number of particles generated at each observational moments
#(M: total number of observational time points)
#(N: total number of vertices)
#thr1: stopping threshold for EM iteration
#thr2: stopping threshold for Robbins Monro iteration
#process: "ER" or "PR"
#init: initialized value for p, q, gamma, alpha, beta

#step_init: 
#n1: number of sample paths used to estimate the derivative of expectation
#n2: the number of updates of parameters before averaging
#times:

#MaxIter1: maximum number of iterations allowed for EM 
#MaxIter2: 
#prop_max: accept proposed sample path when prop_max exceeded
#clusterSizePF: number of cores for particle filtering (exact division of M)
#clusterSizeRM: number of cores for RM (exact division of M)
#particleParallel: TRUE if running particle filtering in parallel
#outFile: name of outFile
#------Output---------------------------------
##parameter estimation
#p: birth rate
#q: death rate 
#gamma: rate parameter 
#alpha: birth rate 
#beta: death rate

EM_parallel = function(Y, t_obs, B, process = "ER", thr=1e-1, init = NULL, D = 50, MaxIter = 300, 
                       prop_max = 100, burnIn = 50, phi, H = 100, particleParallel = FALSE, 
                       MHParallel = FALSE, outFile = "output.txt"){
  if (length(init)==0){
    init = rep(0.1,5)
  }
  p_cur = init[1]
  q_cur = init[2]
  gamma_cur = init[3]
  alpha_cur = init[4]
  beta_cur = init[5]
  cat("Initialization:", init)
  cat("Initialization:", init, "\n", sep=" ", file = outFile, append=TRUE)

  M = length(t_obs)
  #stopping condition for EM:
  stopConEM = 1
  iter = 0

  #keep track of each EM updates??? ####################################
  while (stopConEM > thr & iter < MaxIter){
    cat("---------------------------\n")
    cat("EM Iteration",iter+1,":Generating particles...\n")
    cat("---------------------------\n", sep = " ", file = outFile, append = TRUE)
    cat("EM Iteration",iter+1,":Generating particles...\n", sep = " ", file = outFile, append = TRUE)
    
    time1 = proc.time()
    if(particleParallel){
      partic = particleFil_parallel (B, Y, t_obs, p_cur, q_cur, gamma_cur, alpha_cur, beta_cur, 
                                     process = "ER", outFile = outFile)
    }else{
      partic = particleFil(B, Y, t_obs, p_cur, q_cur, gamma_cur, alpha_cur, beta_cur, 
                           process = "ER", outFile = outFile)
    }
    
    cat("  Finished particle filtering\n")
    cat("  Finished particle filtering\n", sep = " ", file = outFile, append = TRUE)
    time_used2 = proc.time() - time1
    
    cat("Time elapsed:", time_used2[3],"\n")
    cat("Time elapsed:", time_used2[3],"\n", file = outFile, append = TRUE)

    cat("  Start M-step...\n")
    cat("  Start M-step...\n", sep = " ", file = outFile, append = TRUE)

    # Sample H true network in the hidden layer from B particles
    # then for each of them sample D paths that connects the truth
    # \hat p = E[E[\sum(I{w_{n-1} = 0, w_m = 1)]]/E[E[\sum(I_{w_{n-1} = 0})]]
    counts_vec = rep(0, 5)
    ids = sample(1:B, H)
    cat("   Sampling MH paths for b = ...\n")
    cat("   Sampling MH paths for b = ...\n", sep = " ", file = outFile, append = TRUE)

    for (h in 1:H){
      cat("    ", h, ":", ids[h], "\n")
      cat("    ", h, ":", ids[h], "\n", sep = " ", file = outFile, append = TRUE)
      b = ids[h]
      Y_true = partic$particles[[b]]
      W_true = partic$W[b,]

      # sample D paths that connects the truth, (parallel for M here, not D)
      if (MHParallel){
        clusterSize = as.numeric(Sys.getenv("NSLOTS"))
        if (is.na(clusterSize)) clusterSize = 1

        #initialize a cluster for parallel computing
        cluster = makeCluster(clusterSize, type = "SOCK")

        #sample D paths that connects Y_true, W_true
        results = clusterApply(cluster, 1:(M-1), getSamplePath_parallel, D, Y_true, W_true, t_obs,
                               p_cur, q_cur, gamma_cur, process = "ER", prop_max)
        # counts_vec = c(sum_01, sum_0, sum_10, sum_1, sum_R)
        for (cl in 1:(M-1)){
          counts_vec = counts_vec + results[[cl]]
        }
        stopCluster(cluster)
      }
      else{
        #sample D paths that connects Y_true, W_true
        for (m in 1:(M-1)){
          results = getSamplePath_parallel(m, D, Y_true, W_true, t_obs,
                                           p_cur, q_cur, gamma_cur, process = "ER", prop_max)
          counts_vec = counts_vec + results
        }
      }
    }

    #EM update
    pqgamma = update_pqgamma (counts_vec, H, D, t_obs, M)
    alphaBeta = update_alpha_beta (partic$particles, Y, burnIn, phi)

    p_nex = pqgamma[1]
    q_nex = pqgamma[2]
    gamma_nex = pqgamma[3]
    alpha_nex = alphaBeta[1]
    beta_nex = alphaBeta[2]

    #EM stopping condition
    stopConEM = sqrt(sum((c(p_nex, q_nex, gamma_nex, alpha_nex,
                            beta_nex)-c(p_cur, q_cur, gamma_cur, alpha_cur, beta_cur))^2))/sqrt(sum((c(p_cur, q_cur, gamma_cur, alpha_cur, beta_cur))^2))
    cat("\n Relative change:", stopConEM, "\n")
    cat("\n Relative change:", stopConEM, "\n", sep = " ", file = outFile, append = TRUE)

    p_cur = p_nex
    q_cur = q_nex
    gamma_cur = gamma_nex
    alpha_cur = alpha_nex
    beta_cur = beta_nex
    iter = iter + 1
  }

  cat("---------------------------\n")
  cat("EM Iterations Done!\n")
  cat("---------------------------\n", sep = " ", file = outFile, append = TRUE)
  cat("EM Iterations Done!\n", sep = " ", file = outFile, append = TRUE)

  return (list(Estimation = c(p_cur,q_cur,gamma_cur,alpha_cur, beta_cur),
               iteration = iter))
}

# Y = Y[1:2]
# t_obs = t_obs[1:2]

