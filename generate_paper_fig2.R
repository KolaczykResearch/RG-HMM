source("Model/utilities.R", chdir=T)
source("Model/GetSimSeq.R", chdir=T)

# generate noisy sequence of networks from process "ER" or "PR"
getDataSet = function(N, M, param_val, rate_obs, process = "ER"){
  p = param_val[1]
  q = param_val[2]
  gamma = param_val[3]
  alpha = param_val[4]
  beta = param_val[5]
  t_obs = cumsum(rep(1/rate_obs, M))
  # initialize the starting network with a random graph with only 1 edge
  edges = sample_modified(get_triangular(1:N^2, N), 1)
  net_start = matrix(0, N, N)
  for (i in 1:length(edges)){
    pos = get_pos(edges[i], N)
    net_start[pos[1], pos[2]] = 1
    net_start[pos[2], pos[1]] = 1
  }
  # generate sequence of networks from process "ER" or "PR", add noise to all
  res = get_SimSeq_all(p, q, gamma, N, t_obs, net_start = net_start, process = process)
  net_hidden_noise = add_Noise(res$network_seq_hidden, alpha, beta)
  return(list(t_obs = t_obs, t_all = res$t_transition,
              net_hidden = res$network_seq_hidden, 
              bin_hidden = res$bin_var_hidden,
              net_hidden_noise = net_hidden_noise))
}

set.seed(5963)
N=100; M=220; p=0.9; q=0.1; gamma=2; alpha=0.01; beta=0.01; rate_obs=1.5
data_ER = getDataSet(N, M, c(p, q, gamma, alpha, beta), rate_obs, process='ER')
data_PR = getDataSet(N, M, c(p, q, gamma, alpha, beta), rate_obs, process='PR')

plot((1:length(data_ER$t_all))/N, sapply(data_ER$net_hidden, getSizeGCC, order = 1)/N, type='l',
     xlab='Time steps/Number of vertices', ylab='Proportion of Vertices in Giant Component', 
     col='blue', lty=3, lwd=2)
lines((1:length(data_ER$t_all))/N, sapply(data_ER$net_hidden_noise, getSizeGCC, order = 1)/N, 
      col='green', lty=1, lwd=2)
lines((1:length(data_PR$t_all))/N, sapply(data_PR$net_hidden, getSizeGCC, order = 1)/N,
      col='red', lty=5, lwd=2)
lines((1:length(data_PR$t_all))/N, sapply(data_PR$net_hidden_noise, getSizeGCC, order = 1)/N,
      col='purple', lty=4, lwd=2)
legend('bottomright', 
       c('ER Percolation', 'PR Percolation', 'ER with noise', 'PR with noise'), 
       lty=c(3,5,1,4),
       col=c('blue', 'red', 'green', 'purple'), lwd=2)
# 8.45*6.11

