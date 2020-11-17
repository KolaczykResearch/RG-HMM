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


# # getBayesFactor ()-----------------------------------
# # -------Input----------------------------------
# # Y: Observed sequence of networks (with noise)
# # t_obs: A vetor of observational time points
# # para_est: estimated parameters
# # B: Number of particles generated at each observational moments
# # (M: Total number of observational time points)
# # (N: Total number of vertices)
# # process: "ER" or "PR"
# # D: The length of the chain for mcmc
# # prop_max: accept proposed sample path when prop_max exceeded
# # phi: Number of particles used to compute alpha and gamma
# # H: Number of true network sequences used in mcmc
# # particleParallel: TRUE if running particle filtering in parallel
# # burnIn_mcmc: Number of observational moments dropped in mcmc
# # margParallel: TRUE if calculating marglik in parallel
# # ------Output---------------------------------
# # A number
# getBayesFactor = function(Y, t_obs, para_est_ER, para_est_PR, B, process, D = 10, prop_max = 30, H,
#                           particleParallel = FALSE, burnIn_mcmc = 0, margParallel = T, seed_parallel = NULL){
#   cat('==================================\n')
#   cat('==== Calculating loglik for ER...\n')
#   cat('==================================\n')
#   time1 = proc.time()
#   l_ER = getMargLik2(N, Y, t_obs, para_est_ER, B, process = 'ER', parallel = margParallel)
#   cat(l_ER, '\n')
#   time_used = proc.time() - time1
#   cat("  Time elapsed:", time_used[3],"\n")
# 
#   cat('==================================\n')
#   cat('==== Calculating loglik for PR...\n')
#   cat('==================================\n')
#   time1 = proc.time()
#   l_PR = getMargLik2(N, Y, t_obs, para_est_PR, B, process = 'PR', parallel = margParallel)
#   cat(l_PR, '\n')
#   time_used = proc.time() - time1
#   cat("  Time elapsed:", time_used[3],"\n")
# 
#   cat('The loglik difference is:', l_ER-l_PR, '\n')
# 
#   cat('==================================\n')
#   cat('==== Calculating FIM for ER...\n')
#   cat('==================================\n')
#   time1 = proc.time()
#   FIM_ER = getFIM(Y, t_obs, para_est_ER, B, process="ER", D, prop_max, H, particleParallel, burnIn_mcmc)
#   cat('det:', det(FIM_ER), '\n')
#   cat('logdet:', log(det(FIM_ER)), '\n')
#   time_used = proc.time() - time1
#   cat("  Time elapsed:", time_used[3],"\n")
# 
#   cat('==================================\n')
#   cat('==== Calculating FIM for PR...\n')
#   cat('==================================\n')
#   time1 = proc.time()
#   FIM_PR = getFIM(Y, t_obs, para_est_PR, B, process="PR", D, prop_max, H, particleParallel, burnIn_mcmc)
#   cat('det:', det(FIM_PR), '\n')
#   cat('logdet:', log(det(FIM_PR)), '\n')
#   time_used = proc.time() - time1
#   cat("  Time elapsed:", time_used[3],"\n")
# 
#   cat('\n==== Calculating Bayes Factor...\n')
#   bf = l_ER - l_PR + (-1/2)*(log(det(FIM_ER))-log(det(FIM_PR)))
#   #bf = exp(bf)
#   return(bf)
# }
# 
# # getFIM ()-----------------------------------
# # -------Input----------------------------------
# # Y: Observed sequence of networks (with noise)
# # t_obs: A vetor of observational time points
# # para_est: estimated parameters
# # B: Number of particles generated at each observational moments
# # (M: Total number of observational time points)
# # (N: Total number of vertices)
# # process: "ER" or "PR"
# # D: The length of the chain for mcmc
# # prop_max: accept proposed sample path when prop_max exceeded
# # phi: Number of particles used to compute alpha and gamma
# # H: Number of true network sequences used in mcmc
# # particleParallel: TRUE if running particle filtering in parallel
# # burnIn_mcmc: Number of observational moments dropped in mcmc
# # ------Output---------------------------------
# # 5X5 matrix 
# getFIM = function(Y, t_obs, para_est, B, process, D, prop_max = 50, H, 
#                   particleParallel = FALSE, burnIn_mcmc = 0){
#   M = length(t_obs)
#   N = ncol(Y[[1]])
#   edgeIndex = matrix(1:(N^2), N, N)[upper.tri(matrix(1:(N^2), N, N))]
#   
#   # generate B particles 
#   cat("--> Generating particles...\n")
#   time1 = proc.time()
#   if(particleParallel){
#     partic = particleFil_parallel (B, Y, t_obs, para_est[1], para_est[2], para_est[3], para_est[4], para_est[5],
#                                    process = process, edgeIndex)
#   }else{
#     partic = particleFil(B, Y, t_obs, para_est[1], para_est[2], para_est[3], para_est[4], para_est[5],
#                          process = process, E_step = F, edgeIndex)
#   }
#   cat("  Finished particle filtering. \n")
#   time_used2 = proc.time() - time1
#   cat("  Time elapsed:", time_used2[3],"\n")
#   
#   # calculate derivatives 
#   set = sample(1:B, H, replace = T, prob = partic$weights)
#   # each row: c(sum_a, sum_b, sum_c, sum_d)
#   sum_vec = get_stats_from_partic(partic, Y, set, H, N) 
#   # each row: c(sum_01, sum_0, sum_10, sum_1, sum_R)
#   counts_vec = get_stats_from_mcmc(partic, Y, obs, set, H, D, para_est, process, prop_max, burnIn_mcmc)
#   
#   # for numerical instability (to avoild NA), make parameters not equal to zero by adding tiny amount
#   para_est[which(para_est==0)] = 1e-7
#   
#   Hessian = get_expected_Hessian(sum_vec, counts_vec, para_est, D)
#   gradient = get_gradient(sum_vec, counts_vec, para_est, D, H, t_obs, M)
#   gradient_expected = apply(gradient, 2, mean)
#   gradient_prod_expected = get_gradient_prod_expected(gradient)
#   
#   FIM = -Hessian - gradient_prod_expected + outer(gradient_expected, gradient_expected, '*')
#   return(FIM)
# }
# 
# # get sufficient stats for calculating derivatives from particles
# # i.e. phi X 4 matrix, each row: sum_m a_m, sum_m b_m, ... 
# get_stats_from_partic = function(partic, Y, set, H, N){
#   B = length(partic$particles)
#   M = length(partic$particles[[1]])
#   mat = matrix(0, H, 4) 
#   for (id in 1:H){
#     i = set[id]
#     sum_a = 0
#     sum_b = 0
#     sum_c = 0
#     sum_d = 0
#     temp = get_particle_path(partic, i, M)
#     Y_true = temp$Y_true
#     Y_true = lapply(Y_true, format_back, N=N)
#     for (j in 1:M){
#       sum_a = sum_a + (sum(Y[[j]]==1 & Y_true[[j]]==1)-N)/2
#       sum_b = sum_b + sum((Y[[j]] - Y_true[[j]])==-1)/2
#       sum_c = sum_c + sum((Y[[j]] - Y_true[[j]])==1)/2
#       sum_d = sum_d + (sum(Y[[j]]==0 & Y_true[[j]]==0)-N)/2
#     }
#     mat[id, ] = c(sum_a, sum_b, sum_c, sum_d)
#   }
#   return(mat)
# }
# 
# # get sufficient stats from augmented paths by mcmc 
# # return HX5 matrix, each row c(sum_01, sum_0, sum_10, sum_1, sum_R)
# get_stats_from_mcmc = function(partic, Y, obs, set, H, D, para_est, process, prop_max, burnIn_mcmc){
#   cat("--> augmenting by MCMC ...\n")
#   time1 = proc.time()
#   # Sample H true network in the hidden layer from B particles
#   # then for each of them sample D paths that connects the truth
#   # list_sum = list()
#   counts_vec = matrix(0, H, 5) # each row: c(sum_01, sum_0, sum_10, sum_1, sum_R)
#   last_particle_index = set
#   for (h in 1:H){
#     index = last_particle_index[h]
#     # for checking purose
#     # list_sum[[h]] = matrix(0, M-1, 5)
#     #if (h%%100 == 0) 
#     cat(h, ", ")
#     temp = get_particle_path(partic, index, M)
#     Y_true = temp$Y_true
#     W_true = temp$W_true
#     Y_true = lapply(Y_true, format_back, N=N)
#     # sample D paths that connects Y_true, W_true (could parallel)
#     for (m in 1:(M-1)){
#       results = getSamplePath_nonparallel(m, D, Y_true, W_true, t_obs,
#                                           para_est[1], para_est[2], para_est[3], process = process, prop_max, burnIn_mcmc)
#       counts_vec[h, ] = counts_vec[h, ] + results
#       # list_sum[[h]][m-burnIn, ] = results
#     }
#   }
#   time_used2 = proc.time() - time1
#   cat("  Finished MCMC\n")
#   cat("  Time elapsed:", time_used2[3],"\n")
#   return(counts_vec)
# }
# 
# # return 5X5 matrix with each row for p, q, gamma, alpha, beta
# get_expected_Hessian = function(sum_vec, counts_vec, para_est, D){
#   Hessian = matrix(0, 5, 5)
#   p = para_est[1]
#   q = para_est[2]
#   gamma = para_est[3]
#   alpha = para_est[4]
#   beta = para_est[5]
#   Hessian[1, 1] = (-1/p^2)*mean(counts_vec[ ,1]/D) - (1/(1-p)^2)*mean(counts_vec[ ,2]/D)
#   Hessian[2, 2] = (-1/q^2)*mean(counts_vec[ ,3]/D) - (1/(1-q)^2)*mean(counts_vec[ ,4]/D)
#   Hessian[3, 3] = (-1/gamma^2)*mean(counts_vec[ ,5]/D)
#   Hessian[4, 4] = (-1/alpha^2)*mean(sum_vec[ ,3]) - (1/(1-alpha)^2)*mean(sum_vec[ ,4])
#   Hessian[5, 5] = (-1/beta^2)*mean(sum_vec[ ,2]) - (1/(1-beta)^2)*mean(sum_vec[ ,1])
#   return(Hessian)
# }
# 
# # return phi X 5 matrix, each row is gradient given one sampled truth
# get_gradient = function(sum_vec, counts_vec, para_est, D, H, t_obs, M){
#   gradient = matrix(0, H, 5) # each row is the gradient for theta given one sampled truth
#   p = para_est[1]
#   q = para_est[2]
#   gamma = para_est[3]
#   alpha = para_est[4]
#   beta = para_est[5]
#   gradient[, 1] = (1/p)*counts_vec[ ,1]/D - (1/(1-p))*counts_vec[ ,2]/D
#   gradient[, 2] = (1/q)*counts_vec[ ,3]/D - (1/(1-q))*counts_vec[ ,4]/D
#   gradient[, 3] = -(t_obs[M]-t_obs[1]) + counts_vec[ ,5]/D
#   gradient[, 4] = (1/alpha)*sum_vec[ ,3] - (1/(1-alpha))*sum_vec[ ,4]
#   gradient[, 5] = (1/beta)*sum_vec[ ,2] - (1/(1-beta))*sum_vec[ ,1]
#   return(gradient)
# }
# 
# get_gradient_prod_expected = function(gradient){
#   res = matrix(0, 5, 5)
#   for(i in 1:5){
#     for (j in 1:5){
#       res[i, j] = mean(gradient[, i]*gradient[, j])
#     }
#   }
#   return(res)
# }
# 
# 
# 
# 
# 
