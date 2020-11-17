source("Model/DataSet.R", chdir=T)
library(ggplot2)

# simulate n_lines sequence according to ER/PR, return size,mean,se of 1st/2nd cc at each obs time points
get_sequences = function(M, N, p, q, gamma, alpha, beta, rate_obs, process, n_lines){
  param_val = c(p, q, gamma, alpha, beta)
  order1 = matrix(0, n_lines, M)
  # order2 = matrix(0, n_lines, M)
  for(i in 1:n_lines){
    data = DataSet (N, M, param_val, rate_obs, process)
    t_obs = data$t_obs
    order1[i, ] = sapply(data$net_obs_noise, getSizeGCC, order = 1)/N
    # order2[i, ] = sapply(data$net_obs_noise, getSizeGCC, order = 2)
    if (i%%50==0) cat(i, ',')
  }
  #return(list('order1'=order1, 'order2'=order2))
  return(order1)
}

# get data frame for ggplot
get_df = function(M, N, p, q, gamma, alpha, beta, rate_obs, n_lines){
  t_obs = cumsum(rep(1/rate_obs, M))
  ER = get_sequences(M, N, p, q, gamma, alpha, beta, rate_obs, process='ER', n_lines)
  PR = get_sequences(M, N, p, q, gamma, alpha, beta, rate_obs, process='PR', n_lines)
  mean = c(apply(ER, 2, mean), apply(PR, 2, mean))
  sd = c(apply(ER, 2, sd), apply(PR, 2, sd))
  process = rep(c('ER', 'PR'), each=M)
  t_obs = rep(t_obs, 2)
  N_v = rep(N, 2*M)
  df = data.frame(process, t_obs, mean, sd, N_v)
  return(df)
}

get_df_scaled = function(N, x_lim=1.5, p, q, gamma, alpha, beta, rate_obs, n_lines){
  # the # of obs. time points needed for scaled_time to reach 2 with rate_obs
  # M = scaled_time_i*N*rate_obs
  M = ceiling(x_lim*N*rate_obs)
  t_obs = cumsum(rep(1/rate_obs, M))
  scaled_t_obs = t_obs/N*gamma
  ER = get_sequences(M, N, p, q, gamma, alpha, beta, rate_obs, process='ER', n_lines)
  PR = get_sequences(M, N, p, q, gamma, alpha, beta, rate_obs, process='PR', n_lines)
  mean = c(apply(ER, 2, mean), apply(PR, 2, mean))
  sd = c(apply(ER, 2, sd), apply(PR, 2, sd))
  process = rep(c('ER', 'PR'), each=M)
  # t_obs = rep(t_obs, 2)
  scaled_t_obs = rep(scaled_t_obs, 2)
  N_v = rep(N, 2*M)
  df = data.frame(process, scaled_t_obs, mean, sd, N_v)
  return(df)
}

RATE_VAL = 1.5
set.seed(12)
N=20; p=0.9; q=0.1; gamma=2; alpha=0.01; beta=0.01
rate_obs = RATE_VAL
df10 = get_df_scaled(N, x_lim=1.5, p, q, gamma, alpha, beta, rate_obs, n_lines=500)
plt10 = ggplot(df10, aes(x=scaled_t_obs, y=mean, group=process, color=process, shape=process)) +
  geom_line() +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  ggtitle(bquote(atop(~p==.(p)~','~q==.(q)~','~gamma==.(gamma)~','~
                        alpha==.(alpha)~','~beta==.(beta)~','~kappa==.(rate_obs)))) +
  labs(y='Proportion of vertices in GCC', x='Normalized Time') +
  geom_vline(xintercept=c(0.8, 1.34, 1.86, 2.4), size=1, color='darkgray') 
plt10 = plt10 + scale_x_continuous(breaks=c(0.8, 1.34, 1.86, 2.4)) +
                scale_shape_manual(values=c(1,2)) + theme_bw()
plt10
# 6.89*2.94
