source("utilities.R")
source("prob_fun.R")

#return the shortest path that brings y1,w1 to y2,w2
#w1 is actually redundant input 
#return netSeq_path, binSeq_path, with the first element y1,w1, last element y2,w2
get_shortestPath = function(y1, w1, y2, w2){
  N = ncol(y1)
  diffMat = y2-y1
  zeroToOne = get_triangular(which(diffMat == 1), N) #0->1
  oneToZero = get_triangular(which(diffMat == -1),N) #1->0
  #cat("Sampling path...")
  
  ##IF y1, y2 are the same
  if (length(zeroToOne) == 0 & length(oneToZero) ==0){
    #cat("ERROR1==")

if (is.empty(y2)){
  w2 = 0
}
if (is.full(y2)){
  w2 = 1
}
    netSeq_path = list()
    binSeq_path = c(w1,1-w2,w2)
    netSeq_path[[1]] = y1
    netSeq_path[[2]] = y1
    netSeq_path[[3]] = y2
    if(w2==0){
      pos = get_pos(sample_modified(c(get_triangular(which(y1 == 0), N)),1), N)
      netSeq_path[[2]][pos[1],pos[2]] = 1
      netSeq_path[[2]][pos[2],pos[1]] = 1
    }else{
      pos = get_pos(sample_modified(c(get_triangular(which(y1 == 1), N)),1), N)
      netSeq_path[[2]][pos[1],pos[2]] = 0
      netSeq_path[[2]][pos[2],pos[1]] = 0
    }
  }
  
  #Case 2
  if (length(zeroToOne) == 0 & length(oneToZero) != 0){
    #cat("ERROR2==")
    netSeq_path = list()
#############
if (w2==1 & is.empty(y2)){
  w2 = 0
}
#############
    if (w2 == 1){
      
      binSeq_path = rep(c(w1,0,0,1),c(1,length(oneToZero),1,1))
      netSeq_path[[1]] = y1
      index = 1
      for (i in oneToZero){
        pos = get_pos(i, N)
        netSeq_path[[index+1]] = netSeq_path[[index]]
        netSeq_path[[index+1]][pos[1],pos[2]] = 0
        netSeq_path[[index+1]][pos[2],pos[1]] = 0
        index = index+1
      }
      #randomly delete one edge then add
      pos = get_pos(sample_modified(c(get_triangular(which(netSeq_path[[index]] == 1), N)),1), N)
      netSeq_path[[index + 1]] = netSeq_path[[index]]
      netSeq_path[[index+1]][pos[1],pos[2]] = 0
      netSeq_path[[index+1]][pos[2],pos[1]] = 0
      index = index+1
      netSeq_path[[index + 1]] = netSeq_path[[index-1]] 
      
    }else{
      
      binSeq_path = rep(c(w1,0),c(1,length(oneToZero)))
      netSeq_path[[1]] = y1
      index = 1
      for (i in oneToZero){
        pos = get_pos(i, N)
        netSeq_path[[index+1]] = netSeq_path[[index]]
        netSeq_path[[index+1]][pos[1],pos[2]] = 0
        netSeq_path[[index+1]][pos[2],pos[1]] = 0
        index = index+1
      }
    }
  }
  
  #Case 3
  if(length(oneToZero) == 0 & length(zeroToOne) != 0){
    #cat("ERROR3==")
    netSeq_path = list()
#############
if (w2==0 & is.full(y2)){
  w2 = 1
}
#############
    if (w2 == 0){
      
      binSeq_path = rep(c(w1,1,1,0),c(1,length(zeroToOne),1,1))
      netSeq_path[[1]] = y1
      index = 1
      for (i in zeroToOne){
        pos = get_pos(i, N)
        netSeq_path[[index+1]] = netSeq_path[[index]]
        netSeq_path[[index+1]][pos[1],pos[2]] = 1
        netSeq_path[[index+1]][pos[2],pos[1]] = 1
        index = index+1
      }
      #randomly add one edge then delete
      pos = get_pos(sample_modified(c(get_triangular(which(netSeq_path[[index]] == 0), N)),1), N)
      netSeq_path[[index + 1]] = netSeq_path[[index]]
      netSeq_path[[index+1]][pos[1],pos[2]] = 1
      netSeq_path[[index+1]][pos[2],pos[1]] = 1
      index = index+1
      netSeq_path[[index + 1]] = netSeq_path[[index-1]] 
      
    }else{
      
      binSeq_path = rep(c(w1,1),c(1,length(zeroToOne)))
      netSeq_path[[1]] = y1
      index = 1
      for (i in zeroToOne){
        pos = get_pos(i, N)
        netSeq_path[[index+1]] = netSeq_path[[index]]
        netSeq_path[[index+1]][pos[1],pos[2]] = 1
        netSeq_path[[index+1]][pos[2],pos[1]] = 1
        index = index+1
      }
    }
  }
  
  #Case 4
  if(length(oneToZero) != 0 & length(zeroToOne) != 0){
    #cat("ERROR4==")
    netSeq_path = list()
    if (w2 == 0){
      
      binSeq_path = rep(c(w1,1,0),c(1,length(zeroToOne),length(oneToZero)))
      netSeq_path[[1]] = y1
      index = 1
      for (i in zeroToOne){
        pos = get_pos(i, N)
        netSeq_path[[index+1]] = netSeq_path[[index]]
        netSeq_path[[index+1]][pos[1],pos[2]] = 1
        netSeq_path[[index+1]][pos[2],pos[1]] = 1
        index = index+1
      }
      for (i in oneToZero){
        pos = get_pos(i, N)
        netSeq_path[[index+1]] = netSeq_path[[index]]
        netSeq_path[[index+1]][pos[1],pos[2]] = 0
        netSeq_path[[index+1]][pos[2],pos[1]] = 0
        index = index+1
      }
      
    }else{
      
      binSeq_path = rep(c(w1, 0, 1),c(1, length(oneToZero), length(zeroToOne)))
      netSeq_path[[1]] = y1
      index = 1
      for (i in oneToZero){
        pos = get_pos(i, N)
        netSeq_path[[index+1]] = netSeq_path[[index]]
        netSeq_path[[index+1]][pos[1],pos[2]] = 0
        netSeq_path[[index+1]][pos[2],pos[1]] = 0
        index = index+1
      }
      for (i in zeroToOne){
        pos = get_pos(i, N)
        netSeq_path[[index+1]] = netSeq_path[[index]]
        netSeq_path[[index+1]][pos[1],pos[2]] = 1
        netSeq_path[[index+1]][pos[2],pos[1]] = 1
        index = index+1
      }
    }
  }
  if (length(netSeq_path)==2){
    netSeq_path[[4]] = netSeq_path[[2]]
    netSeq_path[[3]] = netSeq_path[[1]]
    netSeq_path[[2]] = netSeq_path[[1]]
    if (is.empty(netSeq_path[[1]])){
      nonEdge = get_triangular(which(netSeq_path[[1]] == 0),N)
      cand = sample_modified(nonEdge,1)
      pos = get_pos(cand,N)
      netSeq_path[[2]][pos[1],pos[2]] = 1
      netSeq_path[[2]][pos[2],pos[1]] = 1
    }else{
      edge = get_triangular(which(netSeq_path[[1]] == 1),N)
      cand = sample_modified(edge,1)
      pos = get_pos(cand,N)
      netSeq_path[[2]][pos[1],pos[2]] = 0
      netSeq_path[[2]][pos[2],pos[1]] = 0
    }
    binSeq_path = isPath_GetBinaryVar(netSeq_path)$binSeq_path
  }
  return (list(netSeq_path = netSeq_path, binSeq_path = binSeq_path))
}

#get changed edge indices at each transitioning time from netSeq_path
get_edgeID_change = function(netSeq_path, N, length){
  id = rep(0, length)
  for(i in 2:length){
    id[i] = get_triangular(which(netSeq_path[[i]]-netSeq_path[[i-1]] !=0),N)
  }
  return (id)
}

#update edgeID_change by inserting cand before r1 and after r2
update_edgeID_insert = function(cand, r1, r2, edgeID_change, length){
  edgeID_change_new = rep(0, length+2)
  edgeID_change_new[1:(r1-1)] = edgeID_change[1:(r1-1)]
  edgeID_change_new[r1] = cand 
  edgeID_change_new[(r1+1):(r2+1)] = edgeID_change[r1:r2]
  edgeID_change_new[r2+2] = cand
  edgeID_change_new[(r2+3):(length+2)] = edgeID_change[(r2+1):length]
  return(edgeID_change_new)
}

#get netSeq_path_prop, binSeq_path_prop from 
#proposed edge change sequence edgeID_change_prop, and netSeq_path[[1]]
getPathFrom_EdgeID = function(edgeID_change_prop, netSeq_path, N){
  netSeq_path_prop = list()
  binSeq_path_prop = rep(0,length(edgeID_change_prop))
  netSeq_path_prop[[1]] = netSeq_path[[1]]
  for(i in 2:length(edgeID_change_prop)){
    edge = get_triangular(which(netSeq_path_prop[[i-1]] == 1),N)
    #nonEdge = get_triangular(which(netSeq_path[[i-1]] == 0),N)
    netSeq_path_prop[[i]] = netSeq_path_prop[[i-1]]
    pos = get_pos(edgeID_change_prop[i], N)
    if (is.element(edgeID_change_prop[i], edge)){
      netSeq_path_prop[[i]][pos[1],pos[2]] = 0
      netSeq_path_prop[[i]][pos[2],pos[1]] = 0
      binSeq_path_prop[i] = 0
    }else{
      netSeq_path_prop[[i]][pos[1],pos[2]] = 1
      netSeq_path_prop[[i]][pos[2],pos[1]] = 1
      binSeq_path_prop[i] = 1
    }
  }
  return(list(netSeq_path_prop = netSeq_path_prop, binSeq_path_prop = binSeq_path_prop))
}


#get proposal path from the current path: netSeq_path, binSeq_path 
#by 3 types of small changes
#return: netSeq_path_prop, binSeq_path_prop
get_proposalPath2 = function(netSeq_path, binSeq_path){
  N = nrow(netSeq_path[[1]])
  length = length(netSeq_path)
  edgeID_change = get_edgeID_change(netSeq_path, N, length)
  
  action = sample(c(1,2,3),1)
  #netSeq_path_prop = list()
  
  #paired insertion 
  if (action == 1){
    cand = sample(get_triangular(1:N^2, N),1)
    position = sort(sample_modified(2:(length-1),2)) 
    r1 = position[1]
    r2 = position[2]
    edgeID_change_prop = update_edgeID_insert(cand, r1, r2, edgeID_change, length)
    res = getPathFrom_EdgeID(edgeID_change_prop, netSeq_path, N)
    netSeq_path_prop = res$netSeq_path_prop
    binSeq_path_prop = res$binSeq_path_prop
    #isPath_GetBinaryVar(netSeq_path_prop)
    #cat("Action 1 ")
    if(!isPath_GetBinaryVar(netSeq_path_prop)$isPath){
      stop("Get wrong path by minor change!!!")
    }
  }
  
  #paired deletion 
  if (action == 2){
    edgeID_change_part = edgeID_change[-length(edgeID_change)]
    if (length(edgeID_change_part) == length(unique(edgeID_change_part)) | length(edgeID_change)<=4){
      #no pair to delete, insert instead
      cand = sample(get_triangular(1:N^2, N),1)
      position = sort(sample_modified(2:(length-1),2)) 
      r1 = position[1]
      r2 = position[2]
      edgeID_change_prop = update_edgeID_insert(cand, r1, r2, edgeID_change, length)
      res = getPathFrom_EdgeID(edgeID_change_prop, netSeq_path, N)
      netSeq_path_prop = res$netSeq_path_prop
      binSeq_path_prop = res$binSeq_path_prop
      #cat("Action 2-1 ")
      if(!isPath_GetBinaryVar(netSeq_path_prop)$isPath){
        stop("Get wrong path by minor change!!!")
      }
    }else{
      dupEle = unique(edgeID_change_part[duplicated(edgeID_change_part)])
      cand = sample_modified(dupEle, 1)
      #pos = get_pos(cand, N)
      pos_deleted = which(edgeID_change == cand)[c(1,2)]
      # for (i in pos_deleted[1]:pos_deleted[2]){
      #   netSeq_path[[i]][pos[1],pos[2]] = 0
      #   netSeq_path[[i]][pos[2],pos[1]] = 0
      # }
      # netSeq_path_prop = netSeq_path[-pos_deleted]
      # binSeq_path_prop = binSeq_path[-pos_deleted]
      edgeID_change_prop = edgeID_change[-pos_deleted]
      res = getPathFrom_EdgeID(edgeID_change_prop, netSeq_path, N)
      netSeq_path_prop = res$netSeq_path_prop
      binSeq_path_prop = res$binSeq_path_prop
      #cat("Action 2 ")
      if(!isPath_GetBinaryVar(netSeq_path_prop)$isPath){
        stop("Get wrong path by minor change!!!")
      }
    }
  }
  
  #permutation
  if (action == 3){
    permutedChange = sample_modified(edgeID_change[2:(length-1)], length-2)
    edgeID_change_prop = c(edgeID_change[1], permutedChange, edgeID_change[length])
    res = getPathFrom_EdgeID(edgeID_change_prop, netSeq_path, N)
    netSeq_path_prop = res$netSeq_path_prop
    binSeq_path_prop = res$binSeq_path_prop
    #isPath_GetBinaryVar(netSeq_path_prop)
    #cat("Action 3 ")
    if(!isPath_GetBinaryVar(netSeq_path_prop)$isPath){
      stop("Get wrong path by minor change!!!")
    }
  }
  
  return (list(netSeq_path_prop = netSeq_path_prop, binSeq_path_prop = binSeq_path_prop))
}

# ##paired insertion
# #return the inserted netSeq_path_prop, binSeq_path_prop
# after_pairedInsert = function (netSeq_path, binSeq_path, N, length, edgeID_change){
#   binSeq_path_prop = rep(0, length+2)
#   netSeq_path_prop = list()
#   cand = sample(get_triangular(1:N^2, N),1)
#   pos = get_pos(cand,N)
#   #postion = sort(sample(2:(length-1),2))
#   #r1 = postion[1]
#   #r2 = postion[2]
#   r1 = sample_modified(2:(length),1)
# 
#   for (i in 1:(r1-1)){
#     netSeq_path_prop[[i]] = netSeq_path[[i]]
#     binSeq_path_prop[i] = binSeq_path[i]
#   }
#   netSeq_path_prop[[r1]] = netSeq_path[[r1-1]]
#   #decide if we should add or delete cand before r1, if cand is not an edge, then add it
#   #and delete right after; otherwise, delete it then add before r1
#   allEdge = get_triangular(which(netSeq_path[[r1-1]] == 1),N)
#   if (is.element(cand, allEdge)){
# 
#     netSeq_path_prop[[r1]][pos[1],pos[2]] = 0
#     netSeq_path_prop[[r1]][pos[2],pos[1]] = 0
#     netSeq_path_prop[[r1+1]] = netSeq_path[[r1-1]]
#     binSeq_path_prop[r1] = 0
#     binSeq_path_prop[r1+1] = 1
#     for(i in (r1+2):(length+2)){
#       netSeq_path_prop[[i]] = netSeq_path[[i-2]]
#       binSeq_path_prop[i] = binSeq_path[i-2]
#     }
#   }else{
#     netSeq_path_prop[[r1]][pos[1],pos[2]] = 1
#     netSeq_path_prop[[r1]][pos[2],pos[1]] = 1
#     netSeq_path_prop[[r1+1]] = netSeq_path[[r1-1]]
#     binSeq_path_prop[r1] = 1
#     binSeq_path_prop[r1+1] = 0
#     for(i in (r1+2):(length+2)){
#       netSeq_path_prop[[i]] = netSeq_path[[i-2]]
#       binSeq_path_prop[i] = binSeq_path[i-2]
#     }
#   }
#   return (list(netSeq_path_prop = netSeq_path_prop,
#                binSeq_path_prop = binSeq_path_prop))
# }
# 
# ##paired deletion
# #return the deleted netSeq_path_prop, binSeq_path_prop
# after_pairedDelete = function (netSeq_path, binSeq_path, N, length, edgeID_change){
# 
# }
# 
# #get proposal path from the current path: netSeq_path, binSeq_path 
# #by 3 types of small changes
# #return: netSeq_path_prop, binSeq_path_prop
# get_proposalPath = function(netSeq_path, binSeq_path){
#   N = nrow(netSeq_path[[1]])
#   length = length(netSeq_path)
#   edgeID_change = get_edgeID_change(netSeq_path, N, length)
#   
#   #action = sample(c(1,2,3),1)
#   action = sample(c(1,3),1)
#   netSeq_path_prop = list()
#   
#   #paired insertion 
#   if (action == 1){
#     res = after_pairedInsert(netSeq_path, binSeq_path, N, length, edgeID_change)
#     netSeq_path_prop = res$netSeq_path_prop
#     binSeq_path_prop = res$binSeq_path_prop
#     #isPath_GetBinaryVar(netSeq_path_prop)
#     cat("1")
#   }
#   
#   #paired deletion 
#   if (action == 2){
#     if (length(edgeID_change) == length(unique(edgeID_change)) | length(edgeID_change)<=3){
#       res = after_pairedInsert(netSeq_path, binSeq_path, N, length, edgeID_change)
#       netSeq_path_prop = res$netSeq_path_prop
#       binSeq_path_prop = res$binSeq_path_prop
#       cat("2-1")
#     }else{
#       dupEle = unique(edgeID_change[duplicated(edgeID_change)])
#       cand = sample_modified(dupEle, 1)
#       pos = get_pos(cand, N)
#       pos_deleted = which(edgeID_change == cand)[c(1,2)]
#       for (i in pos_deleted[1]:pos_deleted[2]){
#         netSeq_path[[i]][pos[1],pos[2]] = 0
#         netSeq_path[[i]][pos[2],pos[1]] = 0
#       }
#       netSeq_path_prop = netSeq_path[-pos_deleted]
#       binSeq_path_prop = binSeq_path[-pos_deleted]
#       cat("2")
#     }
#   }
#   
#   #permutation
#   if (action == 3){
#     
#     if (length(edgeID_change) == length(unique(edgeID_change))){
#       res = after_pairedInsert(netSeq_path, binSeq_path, N, length, edgeID_change)
#       netSeq_path_prop = res$netSeq_path_prop
#       binSeq_path_prop = res$binSeq_path_prop
#       cat("3-1")
#     }else{
#       dupEle = unique(edgeID_change[duplicated(edgeID_change)])
#       if (length(dupEle) == 1){
#         
#         # #delete such pair
#         # cand = sample_modified(dupEle, 1)
#         # pos = get_pos(cand, N)
#         # pos_deleted = which(edgeID_change == cand)[c(1,2)]
#         # for (i in pos_deleted[1]:pos_deleted[2]){
#         #   netSeq_path[[i]][pos[1],pos[2]] = 0
#         #   netSeq_path[[i]][pos[2],pos[1]] = 0
#         # }
#         # netSeq_path_prop = netSeq_path[-pos_deleted]
#         # binSeq_path_prop = binSeq_path[-pos_deleted]
#         
#         res = after_pairedInsert(netSeq_path, binSeq_path, N, length, edgeID_change)
#         netSeq_path_prop = res$netSeq_path_prop
#         binSeq_path_prop = res$binSeq_path_prop
#         cat("3-2-1")
#         
#       }else{
#         #only permutate duplicated elements
#         #sample(dupEle, length(dupEle), replace = F)
#         ##!!!ADD PERMUTATION HERE
#         res = after_pairedInsert(netSeq_path, binSeq_path, N, length, edgeID_change)
#         netSeq_path_prop = res$netSeq_path_prop
#         binSeq_path_prop = res$binSeq_path_prop
#         cat("3-3-1")
#       }
#     }
#   }
#   
#   return (list(netSeq_path_prop = netSeq_path_prop, binSeq_path_prop = binSeq_path_prop))
# }





#getSamplePath()--------------------------------
#get D sample path (s1,v1) that brings y(t1),w(t1) to y(t2),w(t2) conditional on Y(t1), W(t1)
#by Metropolis Hastings
#-------Input----------------------------------
#D: number of sample paths 
#y1, y2: true network in the hidden layer at obs. moments t1, t2
#w1, w2: binary variables at t1, t2
#t1, t2: two obs. moments
#p: birth rate
#q: death rate 
#gamma: rate parameter 
#process: "ER" or "PR"
#(M: total number of observational time points)
#(N: total number of vertices)
#------Output---------------------------------
##get D sample paths (s_1^{(d)}, v_1^{(d)})_{d=1}^D that brings y1,w1 to y2, w2
getSamplePath = function(D, y1, w1, y2, w2, t1, t2, p, q, gamma, process = "ER", prop_max = 100){
  path_shortest = get_shortestPath(y1, w1, y2, w2)
  netSeq_path_cur = path_shortest$netSeq_path
  binSeq_path_cur = path_shortest$binSeq_path
  prop_num = 1
  path_num = 1
  netSeq_path_D = array(list(),D)
  binSeq_path_D = list()
  while (path_num <= D){
    temp = get_proposalPath2 (netSeq_path_cur, binSeq_path_cur)
    prop_num = prop_num + 1
    netSeq_path_prop = temp$netSeq_path_prop 
    binSeq_path_prop = temp$binSeq_path_prop
    
    if(!isPath_GetBinaryVar(netSeq_path_prop)$isPath){
      stop("Get wrong proposed path!!!")
    }
    accept = min(1,exp(getPathProb(netSeq_path_prop, binSeq_path_prop, t1,t2, p, q, gamma))/
                   exp(getPathProb(netSeq_path_cur, binSeq_path_cur, t1,t2, p, q, gamma)))
###########
if(is.nan(accept)){
  accept = 0
}
##########
    if (sample(c(1,0),1, prob = c(accept, 1-accept)) == 1 | prop_num > prop_max){
      #cat(" Accept! ")
      netSeq_path_D[[path_num]] = netSeq_path_prop
      binSeq_path_D[[path_num]] = binSeq_path_prop
      netSeq_path_cur = netSeq_path_prop
      binSeq_path_cur = binSeq_path_prop
      prop_num = 0
      path_num = path_num + 1
      if(!isPath_GetBinaryVar(netSeq_path_cur)$isPath){
        stop("Get wrong path!!!")
      }
    }
  }
  #cat("Got D Paths!\n")
  return (list(netSeq_path_D = netSeq_path_D, binSeq_path_D = binSeq_path_D))
}


# set.seed(7)
# y1 = net_hidden_all[[10]]
# y2 = net_hidden_all[[15]]
# w1 = checkAdding(net_hidden_all[[9]],net_hidden_all[[10]])
# w2 = checkAdding(net_hidden_all[[14]],net_hidden_all[[15]])
# t1 = t_all[10]
# t2 = t_all[15]
# D = 5
# res = getSamplePath (D, y1, w1, y2, w2, t1, t2, p, q, gamma, process = "ER", prop_max = 20)
# res$netSeq_path_D[[1]] 
# res$binSeq_path_D

# netSeq_path = net_hidden_all[2:7]
# binSeq_path = isPath_GetBinaryVar(netSeq_path)$binSeq_path
# getPathProb(netSeq_path, binSeq_path, t1 = t_obs[2],t2 = t_obs[7], p=0, q=1, gamma=0)

#wrapper function of getSamplePath for parallelizaiton
# getSamplePath_parallel = function(processID, D, Y_true, W_true, t_obs,
#                                   p_cur_RM, q_cur_RM, gamma_cur_RM, process = "ER", prop_max, clusterSize){
#   source("MHSamples.R")
#   M = length(t_obs)
#   portion = M/clusterSize
#   score_m_p_natural = 0
#   score_m_q_natural = 0
#   score_m_gamma_natural = 0
#   score_m_grad_p_natural = 0
#   score_m_grad_q_natural = 0
#   score_m_grad_gamma_natural = 0
#
#   for (m in ((processID-1)*portion+1):(processID*portion)){
#     if (m==1) next
#     MHsamples = getSamplePath(D, y1 = Y_true[[m-1]], w1 = W_true[m-1],
#                               y2 = Y_true[[m]], w2 = W_true[m],
#                               t1 = t_obs[m-1], t2 = t_obs[m], p_cur_RM, q_cur_RM, gamma_cur_RM,
#                               process = "ER", prop_max = prop_max)
#     for (d in 1:D){
#       score_m_natural = get_score_m_natural(MHsamples$netSeq_path_D[[d]], MHsamples$binSeq_path_D[[d]],
#                                             t1 = t_obs[m-1], t2 = t_obs[m], p_cur_RM, q_cur_RM, gamma_cur_RM)
#       score_m_p_natural = score_m_p_natural+score_m_natural[1]
#       score_m_q_natural = score_m_q_natural+score_m_natural[2]
#       score_m_gamma_natural = score_m_gamma_natural+score_m_natural[3]
#       score_m_grad_p_natural = score_m_grad_p_natural+score_m_natural[4]
#       score_m_grad_q_natural = score_m_grad_q_natural+score_m_natural[5]
#       score_m_grad_gamma_natural = score_m_grad_gamma_natural+score_m_natural[6]
#     }
#   }
#   return (list(score_m_natural = c(score_m_p_natural, score_m_q_natural, score_m_gamma_natural),
#                score_m_grad_natural = c(score_m_grad_p_natural, score_m_grad_q_natural, score_m_grad_gamma_natural)))
# }

getSamplePath_parallel = function(processID, D, Y_true, W_true, t_obs, p_cur, q_cur, 
                                  gamma_cur, process = "ER", prop_max){
    source("MHSamples.R")
    M = length(t_obs)
    sum_01 = 0
    sum_0 = 0
    sum_10 = 0
    sum_1 = 0
    sum_R = 0 
    
    m = processID + 1
    MHsamples = getSamplePath(D, y1 = Y_true[[m-1]], w1 = W_true[m-1], y2 = Y_true[[m]], w2 = W_true[m], 
                              t1 = t_obs[m-1], t2 = t_obs[m], p_cur, q_cur, gamma_cur, 
                              process = "ER", prop_max = prop_max)
    for (d in 1:D){
      counts = getCount(MHsamples$binSeq_path_D[[d]])
      sum_01 = sum_01 + counts[1]
      sum_0 = sum_0 + counts[2]
      sum_10 = sum_10 + counts[3]
      sum_1 = sum_1 + counts[4]
      sum_R = sum_R + counts[5]
    }
    return (c(sum_01, sum_0, sum_10, sum_1, sum_R))
}

# get count statistics from a series of MH samples
getCount = function(binSeq_path_D){
  len = length(binSeq_path_D)
  sum_R = len - 2
  sum_1 = sum(binSeq_path_D[1:(len-1)])
  sum_0 = len - 1 - sum_1
  diff = binSeq_path_D[1:(len-1)] - binSeq_path_D[2:len]
  sum_10 = sum(diff == 1)
  sum_01 = sum(diff == (-1))
  return (c(sum_01, sum_0, sum_10, sum_1, sum_R))
}


