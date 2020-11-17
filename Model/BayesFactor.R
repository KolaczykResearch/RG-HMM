# getBayesFactor ()-----------------------------------
# -------Input----------------------------------
# Y: Observed sequence of networks (with noise)
# t_obs: A vetor of observational time points
# para_est_ER: estimated parameters under ER
# para_est_PR: estimated parameters under PR
# B: Number of particles generated at each observational moments
# (M: Total number of observational time points)
# (N: Total number of vertices)
# margParallel: TRUE if calculating marglik in parallel 
# ------Output---------------------------------
# A number 
getBayesFactor = function(N, Y, t_obs, para_est_ER, para_est_PR, B, margParallel = T, seed_parallel=NULL){
  cat('==================================\n')
  cat('==== Calculating loglik for ER...\n')
  cat('==================================\n')
  time1 = proc.time()
  l_ER = getMargLik2(N, Y, t_obs, para_est_ER, B, process = 'ER', parallel = margParallel, seed_parallel)
  cat(l_ER, '\n')
  time_used = proc.time() - time1
  cat("  Time elapsed:", time_used[3],"\n")

  cat('==================================\n')
  cat('==== Calculating loglik for PR...\n')
  cat('==================================\n')
  time1 = proc.time()
  l_PR = getMargLik2(N, Y, t_obs, para_est_PR, B, process = 'PR', parallel = margParallel, seed_parallel)
  cat(l_PR, '\n')
  time_used = proc.time() - time1
  cat("  Time elapsed:", time_used[3],"\n")
  
  cat('The loglik difference is:', l_ER-l_PR, '\n')
  
  cat('\n==== Calculating Bayes Factor...\n')
  bf = l_ER-l_PR
  return(bf)
}

