source("utilities.R")

#get_SimSeq()-----------------------------------
#-------Input----------------------------------
#p: birth rate
#q: death rate 
#gamma: rate parameter 
#N: total number of vertices
#t_obs: a vetor of observational time points
#(M: total number of observational time points)
#------Output---------------------------------
#net_seq_obs: network sequence at observational moments with no noise, i.e. 
#A list of M elements, each containing N tims N adjacency matrix
#network_seq_hidden: network sequence at all transitioning moments
#t_transition
#---------------------------------------
get_SimSeq = function(p, q, gamma, N, t_obs, net_start=NULL, w_start=NULL, process = "ER"){
  M = length(t_obs)
  network_seq_obs = list()
  network_seq_hidden = list()
  #start from an empty network
  if (length(net_start) == 0){
    network_cur = matrix(0, N, N)
  }else{
    network_cur = net_start
  }
  #start from w=1
  if (length(w_start) == 0){
    w = 1
  }else{
    w = w_start
  }
  network_seq_hidden[[1]] = network_cur
  t = 0
  m = 0
  t_transition = 0
  count_hidden = 1
  
  while (m<M){
    
    ###--Get w ~ Bernoulli (p)/Bernoulli (1-q)
    if (is.empty(network_cur)){
      #if graph is empty, p = 1
      w = 1
    }else if (is.full(network_cur)){
      #if graph is full, q = 1
      w = 0
    }else{
      if (w==1){
        w = sample(c(1,0), 1, c(1-q, q), replace = FALSE)
      }else{
        w = sample(c(1,0), 1, c(p, 1-p), replace = FALSE)
      }
    }
    
    ###--Update next transition time with t_incre ~ exp(gamma)
    t_incre = rexp(1, gamma)
    t = t + t_incre
    
    ###--Add edge if w=1, delete edge if w=0
    if (process == "ER"){
      network_nex = network_cur
      if (w == 1){
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
    }
    
    # if (process == "PR"){
    #   network_nex = network_cur
    #   if (w == 1){
    #     cand = get_triangular(which(network_cur == 0), N)
    #     if (length(cand)==1){
    #       pos = get_pos(cand, N)
    #     }else{
    #       edge_cand = sample(cand, 2, replace = FALSE)
    #       edge_cand1 = edge_cand[1]
    #       edge_cand2 = edge_cand[2]
    #       pos1 = get_pos(edge_cand1, N)
    #       pos2 = get_pos(edge_cand2, N)
    #       pos = productRuleCC(pos1, pos2, network_cur, N)
    #     }
    #     network_nex[pos[1],pos[2]] = 1 
    #     network_nex[pos[2],pos[1]] = 1
    #   }else{
    #     cand = get_triangular(which(network_cur == 1), N)
    #     if (length(cand)==1){
    #       pos = get_pos(cand, N)
    #     }else{
    #       pos = get_pos(sample(cand,1), N)
    #     }
    #     network_nex[pos[1],pos[2]] = 0
    #     network_nex[pos[2],pos[1]] = 0
    #   }
    # }

    ##--Update network_seq_hidden
    network_seq_hidden[[count_hidden+1]] = network_nex
    count_hidden = count_hidden + 1
    
    ###--Update network_seq_obs
    #number of time points to be updated as network before transitioning
    num = sum(t_obs<t) 
    #cat (t_obs,"\n",num,"\n")
    if (num > 0){
      # for (i in 1:num){
      #   network_seq_obs[[m+i]] = network_cur
      # }
      network_seq_obs = c(network_seq_obs, lapply(1:num, function(i) network_cur))
      m = m + num
      t_obs = t_obs[-(1:num)]
    }
    
    #check if t_obs[1] needs be updated as network after transitioning when t_obs_cur not empty  
    if (length(t_obs)!=0 & t_obs[1]==t){
      network_seq_obs[[m+1]] = network_nex
      m = m + 1
      t_obs = t_obs[-1]
    }
    
    ###--Update network_cur
    network_cur = network_nex
    t_transition = c(t_transition, t)
  }
  return (list(network_seq_obs = network_seq_obs, network_seq_hidden = network_seq_hidden,
               t_transition=t_transition))
}


#add_Noise()-----------------------------------
#this function adds noise to the network sequence at all observation moments
#-------Input----------------------------------
#net_obs: network seq at observational time points
#alpha: type I error (proportion of edges among non-edges)
#beta: tppe II error (proportion of non-edges among edges)
#------Output---------------------------------
#net_seq_obs_noise: observed network sequence with noise, i.e. 
#A list of M elements, each containing N tims N adjacency matrix
#---------------------------------------
add_Noise = function(net_obs, alpha, beta){
  N = nrow(net_obs[[1]])
  pos =0
  for (i in 1:length(net_obs)){
    #among all nonedges, make alpha% of them edge
    nonedge = get_triangular(which(net_obs[[i]] == 0), N)
    temp = ceiling(length(nonedge)*alpha)
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
    temp = ceiling(length(edge)*beta)
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


# #############################
# ############################
set.seed(7)
#random draw 10 time points uniformly from (0, 20)
t_obs = sort(runif(10,0,20))
#or
#t_obs = c(0.5, 1.3, 4, 4.5, 8, 10)
N = 5
p = 0.7; q = 0.3; gamma = 2; alpha = 0.1; beta = 0.5
res = get_SimSeq (p, q, gamma, N, t_obs)

net_obs = res$network_seq_obs
net_obs_noise = add_Noise(net_obs, alpha, beta)
t_all = res$t_transition
net_hidden_all = res$network_seq_hidden

t_obs # observation moments
net_obs_noise #observed networks with noise

net_obs # noise-free networks at observation moments
t_all #all the transitioning time points in the hidden layer
net_hidden_all #all the true networks at transitioning time points in the hidden layer



