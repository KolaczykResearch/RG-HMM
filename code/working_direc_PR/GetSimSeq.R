# """
# This module defines functions: get_SimSeq_all(), add_Noise(), get_SimSeq()
# """

source("utilities.R")

# This function is called in get_SimSeq_all(). 
# Get next binary variable w_{n+1} given previous w_{n}, y_{n} and parameter p, q
get_w_nex = function(w_cur, network_cur, p, q){
  ###--Get w ~ Bernoulli (p)/Bernoulli (1-q)
  if (is.empty(network_cur)){
    #if graph is empty, p = 1
    w = 1
  }else if (is.full(network_cur)){
    #if graph is full, q = 1
    w = 0
  }else{
    if (w_cur==1){
      w = sample(c(1,0), 1, c(1-q, q), replace = FALSE)
    }else{
      w = sample(c(1,0), 1, c(p, 1-p), replace = FALSE)
    }
  }
  return(w)
}

# This function is called in get_SimSeq_all(). 
# Get network_nex by changing one edge (either add or delete determined by w_nex) 
# in network_cur following process rules
network_change = function(network_cur, w_nex, N, process="ER"){
  if (process == "ER"){
    network_nex = network_cur
    if (w_nex == 1){
      pos = sample_modified(which(network_cur ==0 & upper.tri(network_cur)), 1)
      pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[pos]
      network_nex[pos] = 1 
      network_nex[pos1] = 1
    }else{
      pos = sample_modified(which(network_cur ==1 & upper.tri(network_cur)), 1)
      pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[pos]
      network_nex[pos] = 0
      network_nex[pos1] = 0
    }
  }else if (process == "PR"){
      network_nex = network_cur
      igraph_cur = graph_from_adjacency_matrix(network_cur, mode = "undirected")
      clusters = clusters(igraph_cur)
      if (w_nex == 1){
        cand_two = sample_modified(which(network_cur ==0 & upper.tri(network_cur)), 2)
        vertices = productRuleCC(cand_two[1], cand_two[2], clusters, N, w_nex) 
        network_nex[vertices[1], vertices[2]] = 1 
        network_nex[vertices[2], vertices[1]] = 1 
      }else{
        cand_two = sample_modified(which(network_cur ==1 & upper.tri(network_cur)), 2)
        vertices = productRuleCC(cand_two[1], cand_two[2], clusters, N, w_nex) 
        network_nex[vertices[1], vertices[2]] = 0
        network_nex[vertices[2], vertices[1]] = 0
      }
  }
  return(network_nex)
}

# get_SimSeq_all()-------------------------------
# Get simulated network sequence and binary var sequence, at all transitioning moments. 
# -------Input----------------------------------
# p: Birth rate.
# q: Death rate. 
# gamma: Rate of change parameter.
# N: Total number of vertices.
# t_obs: Vetor of observational time points.
# net_start: The initial network of the process, default is empty graph. 
# w_start: 0 or 1, default is 1.
# process: "ER" or "PR", default it "ER".
# ------Output---------------------------------
# network_seq_obs: Network sequence at observational moments with no noise, i.e. 
#              A list of length(t_obs) elements, each containing N tims N adjacency matrix.
# network_seq_hidden: Network sequence at all transitioning moments, including net_start
# bin_var_obs: Binary variable sequence at observational moments. 
# bin_var_hidden: Binary variable sequence at all transitioning moments, including w0 = 1 by default.
# t_transition: All transitioning moments, including 0 for net_start. 
get_SimSeq_all = function(p, q, gamma, N, t_obs, net_start=NULL, w_start=NULL, process = "ER"){
  M = length(t_obs)
  network_seq_obs = list()
  network_seq_hidden = list()
  bin_var_obs = NULL
  bin_var_hidden = NULL
  # start from an empty network if no initial network provided
  if (length(net_start) == 0){
    network_cur = matrix(0, N, N)
  }else{
    network_cur = net_start
  }
  # start from w_cur=1 if no initial provided
  if (length(w_start) == 0){
    w_cur = 1
  }else{
    w_cur = w_start
  }
  bin_var_hidden = w_cur
  network_seq_hidden[[1]] = network_cur
  t = 0
  m = 0
  t_transition = 0
  count_hidden = 1
  
  while (m<M){
    ###--Get w_nex ~ Bernoulli (p)/Bernoulli (1-q)
    w_nex = get_w_nex(w_cur, network_cur, p, q)
    bin_var_hidden = c(bin_var_hidden, w_nex)
    
    ###--Update next transition time with t_incre ~ exp(gamma)
    t_incre = rexp(1, gamma)
    t = t + t_incre
    t_transition = c(t_transition, t)
    
    ###--Add edge if w_nex=1, delete edge if w_nex=0
    network_nex = network_change(network_cur, w_nex, N, process=process)

    ##--Update network_seq_hidden
    network_seq_hidden[[count_hidden+1]] = network_nex
    count_hidden = count_hidden + 1
    
    ###--Update network_seq_obs
    # number of time points to be updated as network_cur before transitioning
    num = sum(t_obs<t) 
    if (num > 0){
      # for (i in 1:num){
      #   network_seq_obs[[m+i]] = network_cur
      # }
      network_seq_obs = c(network_seq_obs, lapply(1:num, function(i) network_cur))
      bin_var_obs = c(bin_var_obs, rep(w_cur, num))
      m = m + num
      t_obs = t_obs[-(1:num)]
    }
    #check if t_obs[1] needs be updated as network_nex after transitioning when t_obs is not empty  
    if (length(t_obs)!=0 & t_obs[1]==t){
      network_seq_obs[[m+1]] = network_nex
      bin_var_obs = c(bin_var_obs, w_nex)
      m = m + 1
      t_obs = t_obs[-1]
    }
    
    ###--Update network_cur
    network_cur = network_nex
    w_cur = w_nex
  }
  return (list(network_seq_obs = network_seq_obs, network_seq_hidden = network_seq_hidden,
               t_transition=t_transition, 
               bin_var_hidden = bin_var_hidden, bin_var_obs = bin_var_obs))
}


# get_SimSeq()-------------------------------
# Get simulated network sequence and binary var sequence, at only observational moments
# -------Input----------------------------------
# same as get_SimSeq_all()
# ------Output---------------------------------
# network_seq_obs: Network sequence at observational moments with no noise, i.e. 
#              A list of length(t_obs) elements, each containing N tims N adjacency matrix.
# bin_var_obs: Binary variable sequence at observational moments. 
get_SimSeq = function(p, q, gamma, N, t_obs, net_start=NULL, w_start=NULL, process = "ER"){
  M = length(t_obs)
  network_seq_obs = list()
  bin_var_obs = NULL
  # start from an empty network if no initial network provided
  if (length(net_start) == 0){
    network_cur = matrix(0, N, N)
  }else{
    network_cur = net_start
  }
  # start from w_cur=1 if no initial provided
  if (length(w_start) == 0){
    w_cur = 1
  }else{
    w_cur = w_start
  }
  t = 0
  m = 0

  while (m<M){
    ###--Get w_nex ~ Bernoulli (p)/Bernoulli (1-q)
    w_nex = get_w_nex(w_cur, network_cur, p, q)
    
    ###--Update next transition time with t_incre ~ exp(gamma)
    t_incre = rexp(1, gamma)
    t = t + t_incre
    
    ###--Add edge if w_nex=1, delete edge if w_nex=0
    network_nex = network_change(network_cur, w_nex, N, process=process)
    
    ###--Update network_seq_obs
    # number of time points to be updated as network_cur before transitioning
    num = sum(t_obs<t) 
    if (num > 0){
      # for (i in 1:num){
      #   network_seq_obs[[m+i]] = network_cur
      # }
      network_seq_obs = c(network_seq_obs, lapply(1:num, function(i) network_cur))
      bin_var_obs = c(bin_var_obs, rep(w_cur, num))
      m = m + num
      t_obs = t_obs[-(1:num)]
    }
    #check if t_obs[1] needs be updated as network_nex after transitioning when t_obs is not empty  
    if (length(t_obs)!=0 & t_obs[1]==t){
      network_seq_obs[[m+1]] = network_nex
      bin_var_obs = c(bin_var_obs, w_nex)
      m = m + 1
      t_obs = t_obs[-1]
    }
    
    ###--Update network_cur
    network_cur = network_nex
    w_cur = w_nex
  }
  return (list(network_seq_obs = network_seq_obs, bin_var_obs = bin_var_obs))
}


# add_Noise()-----------------------------------
# This function adds noise to the network sequence at all observation moments.
#-------Input----------------------------------
# net_obs: Network seq at observational time points.
# alpha: type I error (proportion of edges among non-edges)
# beta: tppe II error (proportion of non-edges among edges)
#------Output---------------------------------
# net_seq_obs_noise: Observed network sequence with noise, i.e. a list of M elements, 
#                         each containing N tims N adjacency matrix.
#---------------------------------------------
add_Noise = function(net_obs, alpha, beta){
  N = nrow(net_obs[[1]])
  pos =0
  for (i in 1:length(net_obs)){
    #among all nonedges, make alpha% of them edge
    nonedge = get_triangular(which(net_obs[[i]] == 0), N)
    temp = floor(length(nonedge)*alpha)
    #temp = ceiling(length(nonedge)*alpha)
    tochange = sample_modified(nonedge, temp)
    if(length(tochange)!=0){
      for (j in 1:length(tochange)){
        pos = get_pos(tochange[j], N)
        net_obs[[i]][pos[2],pos[1]] = 1
        net_obs[[i]][pos[1],pos[2]] = 1
      }
    }
    #among all edges, make beta% of them nonedges
    edge = get_triangular(which(net_obs[[i]] == 1), N)
    temp = floor(length(edge)*beta)
    #temp = ceiling(length(edge)*beta)
    tochange = sample_modified(edge, temp)
    if(length(tochange)!=0){
      for (j in 1:length(tochange)){
        pos = get_pos(tochange[j], N)
        net_obs[[i]][pos[2],pos[1]] = 0
        net_obs[[i]][pos[1],pos[2]] = 0
      }
    }
  }
  return(net_obs)
}
# 
# add_Noise = function(net_obs, alpha, beta){
#   N = nrow(net_obs[[1]])
#   pos =0
#   for (i in 1:length(net_obs)){
#     #among all nonedges, make alpha% of them edge
#     nonedge = get_triangular(which(net_obs[[i]] == 0), N)
#     bool = sample(c(TRUE, FALSE), length(nonedge), prob = c(alpha, 1-alpha), replace = TRUE)    
#     tochange = nonedge[bool]
#     if(length(tochange)!=0){
#       tochange1 = matrix(1:(N^2), N, N, byrow = TRUE)[tochange]
#       net_obs[[i]][tochange] = 1
#       net_obs[[i]][tochange1] = 1
#     }
#     #among all edges, make beta% of them nonedges
#     edge = get_triangular(which(net_obs[[i]] == 1), N)
#     bool = sample(c(TRUE, FALSE), length(edge), prob = c(beta, 1-beta), replace = TRUE)    
#     tochange = edge[bool]
#     if(length(tochange)!=0){
#       tochange1 = matrix(1:(N^2), N, N, byrow = TRUE)[tochange]
#       net_obs[[i]][tochange] = 0
#       net_obs[[i]][tochange1] = 0
#     }
#   }
#   return(net_obs)
# }
# 
# # ===========================================
# # ===========================================
# set.seed(7)
# #random draw 10 time points uniformly from (0, 20)
# t_obs = sort(runif(10,0,20))
# N = 5
# p = 0.7; q = 0.3; gamma = 2; alpha = 0.1; beta = 0.5
# res = get_SimSeq_all(p, q, gamma, N, t_obs)
# 
# net_obs = res$network_seq_obs
# net_obs_noise = add_Noise(net_obs, alpha, beta)
# bin_obs = res$bin_var_obs
# t_all = res$t_transition
# net_hidden = res$network_seq_hidden
# bin_hidden = res$bin_var_hidden
# 
# t_obs # observation moments
# net_obs_noise #observed networks with noise
# 
# net_obs # noise-free networks at observational moments
# bin_obs # binary variable at observational moments
# t_all # all the transitioning moments in the hidden layer
# net_hidden # true networks at all transitioning moments in the hidden layer
# bin_hidden # binary variable at all transitioning moments
# 
# 
# 



