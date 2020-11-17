#--------utility function----------------------
#---------------------------------------------
#get all indices for the upper triangular elements(diagonal elements excluded)
#(NtimesN matrix) from an array of indices(seq)
# get_triangular = function(seq, N){
#   if(length(seq) ==0){
#     return(NULL)
#   }
#   res = numeric(length(seq))
#   ind = 1
#   for (i in 1:length(seq)){
#     col = ifelse(seq[i]%%N!=0, seq[i]%/%N + 1, seq[i]%/%N)
#     row = ifelse(seq[i]%%N!=0, seq[i]%%N, N)
#     if (col>row){
#       res[ind] = seq[i]
#       ind = ind + 1
#     }
#   }
#   res = res[res!=0]
#   return (res)
# }

get_triangular = function(seq, N){
    M = upper.tri(matrix(0, N, N))
    # return(seq[M[seq]])
    return (seq[which(M[seq])])
}

#covert one-dimensional index to two-dim (row,col) for N times N matrix
get_pos = function (index, N){
    col = ifelse(index%%N!=0, index%/%N + 1, index%/%N)
    row = ifelse(index%%N!=0, index%%N, N)
    return (c(row, col))
}

#check if it is empty graph from adjacency matrix
is.empty = function(mat){
    if (nrow(mat)!=ncol(mat)){
        error("Adjacency matrix must be squared matrix!")
    }
    N = nrow(mat)
    if (sum(mat == 0) == N^2){
        return (TRUE)
    }else{
        return (FALSE)
    }
}

#check if it is full graph from adjacency matrix
is.full = function(mat){
    if (nrow(mat)!=ncol(mat)){
        error("Adjacency matrix must be squared matrix!")
    }
    N = nrow(mat)
    if (sum(mat == 1) == (N^2-N)){
        return (TRUE)
    }else{
        return (FALSE)
    }
}

#rewrite sample function, take n random element from x
sample_modified = function(x, n){
    if (length(x) == 1){
        return(rep(x,n))
    }else{
        return(sample(x,n))
    }
}
#get the number of edges and nonedges of a network
edge_num = function(y1){
    res = sum(y1==1)/2
    return (res)
}

nonedge_num = function(y1){
    N = nrow(y1)
    res = (sum(y1==0)-N)/2
    return (res)
}

#check whether it is adding or deleting edges from net1 to net2
#return 1 if it is adding
# checkAdding = function (net1, net2){
#   #upper.tri() returns the indices
#   uppder_tri_id = upper.tri(net1, diag = F)
#   temp = (net2 - net1)[uppder_tri_id]
#   #check if only one item in temp not equal to zero
#   if (sum(temp!=0)!=1){
#     warning("Infeasible Samples!")
#     return(2)
#   }
#   #add = temp[temp !=0]
#   add = (sum(temp) + 1)/2
#   return (add)
#   # if (add == 1){
#   #   return (1)
#   # }
#   # if (add == (-1)){
#   #   return (0)
#   # }
# }

checkAdding = function (net1, net2){
    net0 = net2 - net1
    if (sum(abs(net0)) > 2){
        warning("Infeasible Samples!")
        return(2)
    }
    add = (sum(net0)/2 + 1)/2
    return (add)
}

#return 1 if it is adding
checkAdding2 = function (net1, net2){
    #upper.tri() returns the indices
    uppder_tri_id = upper.tri(net1, diag = F)
    temp = (net2 - net1)[uppder_tri_id]
    #check if only one item in temp not equal to zero
    if (sum(temp!=0)!=1){
        stop("Infeasible Samples!LALALALA")
        return(2)
    }
    add = temp[temp !=0]
    if (add == 1){
        return (1)
    }
    if (add == (-1)){
        return (0)
    }
}
#get the corresponding binary variables based on
#the network sequence differed by only one edge each time
#and the value(0/1) of binary variable at the first moment
getBinaryVar = function (netSeq, w1){
    len = length(netSeq)
    W = rep(0,len)
    W[1] = w1
    for (i in 2:len){
        W[i] = checkAdding(netSeq[[i-1]], netSeq[[i]])
    }
    return(W)
}

#netSeq = net_hidden_all
#getBinaryVar(netSeq, w1 = 1)

#check if a network sequence is a feasible path from ER process, i.e.
#check if any two consecutive networks differ by only one edge
#AND
##get the corresponding binary variables from the
#feasible network sequence path which differed by only one edge each time
isPath_GetBinaryVar = function (netSeq_path){
    binSeq_path = rep(0, length(netSeq_path))
    binSeq_path[1] = 1
    # for (i in 2:length(netSeq_path)){
    #   res = checkAdding(netSeq_path[[i-1]], netSeq_path[[i]])
    #   #binSeq_path[i] = ifelse(res==1, 1, 0)
    #   binSeq_path[i] = res
    #   if (res == 2){
    #     return(list(binSeq_path = NULL, isPath = FALSE))
    #   }
    # }
    #binSeq_path[i] = ifelse(res==1, 1, 0)
    res = sapply(2:length(netSeq_path), function (i) checkAdding(netSeq_path[[i-1]], netSeq_path[[i]]))
    if (any(res == 2))  return(list(binSeq_path = NULL, isPath = FALSE))
    binSeq_path[2:length(netSeq_path)] = res
    
    return(list(binSeq_path = binSeq_path, isPath = TRUE))
}
#isPath_GetBinaryVar(net_hidden_all)
#isPath_GetBinaryVar(particles[[1]])


#sigmoid
sigmoid = function(x){
    return(1/(1+exp(-x)))
}

sigmoid_dev = function(x){
    return(sigmoid(x)*(1-sigmoid(x)))
}

sigmoid_dev_sec = function(x){
    res = sigmoid(x) - 3*(sigmoid(x))^2 + 2*(sigmoid(x))^3
    return (res)
}

######PR process#################
#get the size of the connected components to which vertex v belongs
getSizeCC = function (network_cur, v, N){
    CC = v
    checkSet = v
    checkedSet = NULL
    while (length(checkSet) !=0){
        neighbors = which(network_cur[checkSet[1],]!=0)
        CC = union(CC, neighbors)
        checkedSet = c(checkedSet, checkSet[1])
        checkSet = checkSet[-1]
        checkSet = setdiff(union(checkSet, neighbors), checkedSet)
    }
    size = length(CC)
    return(size)
}
#network_cur = Y[[1]]
#getSizeCC(network_cur, v=5, N=5)

# determine which edge to add/delete based on product rule of PR process
library(igraph)
productRuleCC = function(pos1, pos2, network_cur, N, w_nex){
  v11 = matrix(rep(1:N, each = N), N, N)[pos1]
  v12 = matrix(rep(1:N, each = N), N, N, byrow = TRUE)[pos1]
  v21 = matrix(rep(1:N, each = N), N, N)[pos2]
  v22 = matrix(rep(1:N, each = N), N, N, byrow = TRUE)[pos2]
  if(w_nex == 1){
    igraph_cur = graph_from_adjacency_matrix(network_cur, mode = "undirected")
    clusters = clusters(igraph_cur)
  }else{
    network_cur[v11, v12] = 0
    network_cur[v12, v11] = 0
    network_cur[v21, v22] = 0
    network_cur[v22, v21] = 0
    igraph_cur = graph_from_adjacency_matrix(network_cur, mode = "undirected")
    clusters = clusters(igraph_cur)
  }
  C11 = clusters$csize[clusters$membership[v11]]
  C12 = clusters$csize[clusters$membership[v12]]
  C21 = clusters$csize[clusters$membership[v21]]
  C22 = clusters$csize[clusters$membership[v22]]
  if (C11*C12 < C21*C22){
    return ( w_nex*c(v11, v12) + (1-w_nex)*c(v21, v22) )
  }else if (C11*C12 == C21*C22){
    temp = sample(c(0,1), 1)
    return ( temp*c(v11, v12) + (1-temp)*c(v21, v22) )
  }else{
    return ( w_nex*c(v21, v22) + (1-w_nex)*c(v11, v12) )
  }
}

# plot the adjacency matrix 
plot_mat = function(y){
  graph = graph_from_adjacency_matrix(y, mode = "undirected")
  plot(graph)
}

# get the nonedge list and the corresponding nonedge product score
get_prod_score_add = function(y, clusters){
  N = nrow(y)
  nonedges = which(y == 0 & upper.tri(y))
  scores = sapply(nonedges, FUN = function(e, N){
    v1 = matrix(rep(1:N, each = N), N, N)[e]
    v2 = matrix(rep(1:N, each = N), N, N, byrow = TRUE)[e]
    C1 = clusters$csize[clusters$membership[v1]]
    C2 = clusters$csize[clusters$membership[v2]]
    return(C1*C2)
  }, N)
  res = cbind(nonedges, scores)
  return(res)
}

# get the edge list and the corresponding edge product score
get_prod_score_delete = function(y, clusters){
  N = nrow(y)
  edges = which(y == 1 & upper.tri(y))
  scores = sapply(edges, FUN = function(e, N){
    v1 = matrix(rep(1:N, each = N), N, N)[e]
    v2 = matrix(rep(1:N, each = N), N, N, byrow = TRUE)[e]
    C1 = clusters$csize[clusters$membership[v1]]
    C2 = clusters$csize[clusters$membership[v2]]
    return(C1*C2)
  }, N)
  res = cbind(edges, scores)
  return(res)
}

# get the size of gcc and 2nd connected component 
getSizeGCC = function(network_cur, order){
  n_ver = ncol(network_cur)
  graph = graph_from_adjacency_matrix(network_cur, mode = "undirected")
  clusters_size = clusters(graph)$csize
  if (order == 1){
    return (max(clusters_size))
  }else if (order == 2){
    n = length(clusters_size)
    if (n > 1){
      id = n-1
      #return (sort(clusters_size)[id] / n_ver)
      return (sort(clusters_size)[id])
    }else{
      return (0)
    }
  }else{
    cat ('Order Not Defined!')
  }
  return(0)
}

# transform data format of network from matrix to array of upper triangle positions 
reformat = function(network, N){
  return(which(network == 1 & upper.tri(matrix(0, N, N))))
}

# transform data format from array to matrix
format_back = function(array, N){
  network = matrix(0, N, N)
  network[array] = 1 
  network = network + t(network)
  return(network)
}

