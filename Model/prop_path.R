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
    cat("Action 1")
  }
  
  #paired deletion 
  if (action == 2){
    edgeID_change_part = edgeID_change[-length(edgeID_change)]
    if (length(edgeID_change_part) == length(unique(edgeID_change_part)) & length(edgeID_change)<=3){
      #no pair to delete, insert instead
      cand = sample(get_triangular(1:N^2, N),1)
      position = sort(sample_modified(2:(length-1),2)) 
      r1 = position[1]
      r2 = position[2]
      edgeID_change_prop = update_edgeID_insert(cand, r1, r2, edgeID_change, length)
      res = getPathFrom_EdgeID(edgeID_change_prop, netSeq_path, N)
      netSeq_path_prop = res$netSeq_path_prop
      binSeq_path_prop = res$binSeq_path_prop
      cat("Action 2-1")
    }else{
      dupEle = unique(edgeID_change_part[duplicated(edgeID_change_part)])
      cand = sample_modified(dupEle, 1)
      pos = get_pos(cand, N)
      pos_deleted = which(edgeID_change == cand)[c(1,2)]
      for (i in pos_deleted[1]:pos_deleted[2]){
        netSeq_path[[i]][pos[1],pos[2]] = 0
        netSeq_path[[i]][pos[2],pos[1]] = 0
      }
      netSeq_path_prop = netSeq_path[-pos_deleted]
      binSeq_path_prop = binSeq_path[-pos_deleted]
      cat("Action 2")
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
    cat("Action 3")
  }
  
  return (list(netSeq_path_prop = netSeq_path_prop, binSeq_path_prop = binSeq_path_prop))
}