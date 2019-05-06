list_sum
partic[[1]][[3008]][1]
partic[[1]][[3008]][2]
partic[[2]][3008, ]
b = 123
test_mcmc(partic[[1]][[b]], partic[[2]][b, ], t_obs, process = "ER", D = 10, prop_max = 30, 
          MHParallel = FALSE, burnIn_mcmc = 1, burnIn = 0, H = 1, p_cur, q_cur, gamma_cur)
test_mcmc(net_obs, partic[[2]][b, ], t_obs, process = "ER", D = 10, prop_max = 30, 
          MHParallel = FALSE, burnIn_mcmc = 1, burnIn = 0, H = 1, p_cur, q_cur, gamma_cur)
test_mcmc(net_obs, bin_obs, t_obs, process = "ER", D = 10, prop_max = 30, 
          MHParallel = FALSE, burnIn_mcmc = 1, burnIn = 0, H = 1, p_cur, q_cur, gamma_cur)
sum(net_obs_noise[[1]] != partic[[1]][[b]][[1]])
sum(net_obs[[1]] != partic[[1]][[b]][[1]])
sum(net_obs[[2]] != partic[[1]][[b]][[2]])
for (m in 1:M){
  cat(sum(net_obs[[m]] != partic[[1]][[b]][[m]])/2, ', ')
}
