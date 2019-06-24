library("MASS")

# help function
vecMatch <- function(x, vec) {
  isTRUE(all.equal(x, vec))
}
# creat the one portion of label_dict from net_obs_noise
# creat target_label_index
# update P1, P21, P3x1
init = function(M, N, net_obs_noise, nState, label_dict, target_label_index, 
                P1, P21, P3x1, m){
  for (i in 2:M){
    network_vec = get_triangular(which(net_obs_noise[[i]] == 1), N)
    # label = graph_to_label(network_vec, N)
    test = sapply(label_dict, vecMatch, vec = network_vec)
    if (any(test)){
      target_label_index[i] = which(test == 1)
    }else{
      nState = nState + 1
      label_dict[[nState]] = network_vec 
      target_label_index[i] = nState
      
    }
    # update P1
    P1[target_label_index[i]] = P1[target_label_index[i]] + 1
    # update P21
    if(i>2){
      P21[target_label_index[i], target_label_index[i-1]] = 
        P21[target_label_index[i], target_label_index[i-1]] + 1
    }
    # update P2x1
    if(i>3){
      if (is.null(P3x1[[target_label_index[i-1]]])){
        P3x1[[target_label_index[i-1]]] = matrix(0, m, m)
      }
      P3x1[[target_label_index[i-1]]][target_label_index[i], target_label_index[i-2]] = 
        P3x1[[target_label_index[i-1]]][target_label_index[i], target_label_index[i-2]] + 1
    }
  }
  return(list(label_dict = label_dict, target_label_index = target_label_index, 
              nState = nState, P1 = P1, P21 = P21, P3x1 = P3x1))
}

# Update label_dict, empirical estimates P1, P21, P3x1 with net_sequence_sampled
update = function(M, N, net_seq_sampled, nState, label_dict, target_label_index, 
                P1, P21, P3x1, m){
  label_index = numeric(M)
  for (i in 2:M){
    network_vec = get_triangular(which(net_seq_sampled[[i]] == 1), N)
    # label = graph_to_label(network_vec, N)
    test = sapply(label_dict, vecMatch, vec = network_vec)
    if (any(test)){
      label_index[i] = which(test == 1)
    }else{
      nState = nState + 1
      label_dict[[nState]] = network_vec 
      label_index[i] = nState
    }
    # update P1
    P1[label_index[i]] = P1[label_index[i]] + 1
    # update P21
    if(i>2){
      P21[label_index[i], label_index[i-1]] = 
        P21[label_index[i], label_index[i-1]] + 1
    }
    # update P2x1
    if(i>3 & label_index[i-1] %in% target_label_index){
      if (is.null(P3x1[[target_label_index[i-1]]])){
        P3x1[[target_label_index[i-1]]] = matrix(0, m, m)
      }
      P3x1[[label_index[i-1]]][label_index[i], label_index[i-2]] = 
        P3x1[[label_index[i-1]]][label_index[i], label_index[i-2]] + 1
    }
  }
  return(list(label_dict = label_dict, nState = nState, P1 = P1, P21 = P21, P3x1 = P3x1))
}

# learnHMM() ----------------------------------
# A spectral algorithm for learning HMM with repeated sampling from hmm 
# with estimated parameters 
# -------Input----------------------------------
# net_obs_noise: observed networks
# (N: Total number of vertices) 
# nChain: the total number of observed sequences to generate from HMM
# lenChain: generating nChain sequence of observed networks with lenChain using para_est
# t_obs: 
# para_est
learnHMM = function(N, net_obs_noise, t_obs, nChain, lenChain = NULL, para_est, process){
  M = length(t_obs)
  label_dict = list()
  target_label_index = numeric(M)
  nState = 0
  chain_fin = 0
  if (is.null(lenChain)){
    m_ = (M-1) * (nChain+1) # observed number of states (upper bound), # of observations
    m = min(m_, 2^choose(N,2))
  }else{
    m_ = M-1 + (lenChain-1)*nChain
    m = min(m_, 2^choose(N,2))
  }
  P1 = numeric(m)
  P21 = matrix(0, m, m)
  P3x1 = array(list(), M)
  P3x1 = lapply(P3x1, FUN = function(x) matrix(0, m, m))
  
  # rename the graph by order of appearance, store it in label_dict
  temp = init(M, N, net_obs_noise, nState, label_dict, target_label_index, P1, P21, P3x1, m)
  label_dict = temp$label_dict
  target_label_index = temp$target_label_index
  nState = temp$nState
  P1 = temp$P1
  P21 = temp$P21
  P3x1 = temp$P3x1
  
  # sample nChain observed sequence from HMM with para_est
  while(chain_fin < nChain){
    cat(' ', chain_fin)
    temp = get_SimSeq(para_est[1], para_est[2], para_est[3], N, t_obs, 
               net_start=net_obs_noise[[1]], w_start=1, process)
    net_seq_sampled = add_Noise(temp$network_seq_obs, para_est[4], para_est[5])
    temp = update(M, N, net_seq_sampled, nState, label_dict, target_label_index, 
                             P1, P21, P3x1, m)
    label_dict = temp$label_dict
    nState = temp$nState
    P1 = temp$P1
    P21 = temp$P21
    P3x1 = temp$P3x1
    chain_fin = chain_fin + 1
  }
  # get P1, P21, P3x1
  P1 = P1[1:nState]/m_
  P21 = P21[1:nState, 1:nState]/(m_-nChain-1)
  P3x1 = lapply(P3x1, FUN = function(x){
    if(!is.null(x)){
      return(x[1:nState, 1:nState]/(m_-2*nChain-2))
    }
  })

  # SVD, compute model parameters
  cat('svd...')
  U = svd(P21)$u
  cat('calculating b1...')
  b1 = t(U)%*%P1
  cat('calculating b_inf...')
  b_inf = ginv(t(P21)%*%U)%*%P1
  cat('calculating b_x...')
  Bx = lapply(P3x1, function(x){
    if(!is.null(x)){
      t(U)%*%x%*%ginv(t(U)%*%P21)
    }else{
      return (NULL)
    }
  })
  return(list(b1 = b1, b_inf = b_inf, Bx = Bx, target_label_index = target_label_index))
}

getMargLik = function(N, net_obs_noise, t_obs, nChain, lenChain = NULL, para_est, process){
  temp = learnHMM(N, net_obs_noise, t_obs, nChain, lenChain, para_est, process)
  M = length(t_obs)
  id_M = temp$target_label_index[M]
  prod = t(temp$b_inf)%*%temp$Bx[[id_M]]
  for(i in (M-1):2){
    id = temp$target_label_index[i]
    prod = prod%*%temp$Bx[[id]]
  }
  prod = prod%*%temp$b1
  return(prod)
}

#unique(temp$target_label_index[-1])

# # map a binary array to decimal
# bin_to_dec = function(x){
#   # x = as.character(as.numeric(x))
#   # b = as.numeric(unlist(strsplit(x, "")))
#   # pow = 2 ^ ((length(b) - 1):0)
#   # return (sum(pow[b == 1]))
#   pow = 2 ^ ((length(x) - 1):0)
#   return (sum(pow[x == 1]))
# } 
# # A function that maps graph to binary array, ie the corresponding label in the state space
# # the graph here is represented as a list of upper triangular indices where edge exists
# graph_to_label = function(network_vec, N){
#   edgeID = matrix(1:N^2, N, N)[upper.tri(matrix(1:N^2, N, N))]
#   binary = rep(0, choose(N,2))
#   # find the corresponding indices in edgeID for network_vec
#   index = match(network_vec, edgeID)
#   binary[index] = 1
#   # label = bin_to_dec(binary)
#   return(binary)
# }

# # get index of observed sequence of networks
# sapply(net_obs_noise, FUN = function(x, label_dict, N){
#   network_vec = get_triangular(which(x == 1), N)
#   state = graph_to_label(network_vec, N)
#   return(match(state, label_dict))
# }, label_dict = label_dict, N = N)

