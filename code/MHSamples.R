source("utilities.R")
source("prob_fun.R")

# Return the shortest path that brings y1,w1 to y2,w2
# w1 is actually redundant input
# Return netSeq_path, binSeq_path, with the first element y1, w1, last element y2, w2
get_shortestPath = function(y1, w1, y2, w2){
    N = ncol(y1)
    diffMat = y2-y1
    # zeroToOne = get_triangular(which(diffMat == 1), N) #0->1
    # oneToZero = get_triangular(which(diffMat == -1),N) #1->0
    zeroToOne = which(diffMat == 1 & upper.tri(diffMat)) #0->1
    oneToZero = which(diffMat == -1 & upper.tri(diffMat)) #1->0
    #cat("Sampling path...")
    
    ##IF y1, y2 are the same
    if (length(zeroToOne) == 0 & length(oneToZero) == 0){
        # if (is.empty(y2)){
        #   w2 = 0
        # }
        # if (is.full(y2)){
        #   w2 = 1
        # }
        netSeq_path = list()
        binSeq_path = c(w1, 1-w2, w2)
        netSeq_path[[1]] = y1
        netSeq_path[[2]] = y1
        netSeq_path[[3]] = y2
        if(w2 == 0){
            pos = sample_modified(which(y1 == 0 & upper.tri(y1)), 1)
            pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[pos]
            netSeq_path[[2]][pos] = 1
            netSeq_path[[2]][pos1] = 1
        }else{
            pos = sample_modified(which(y1 == 1 & upper.tri(y1)), 1)
            pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[pos]
            netSeq_path[[2]][pos] = 0
            netSeq_path[[2]][pos1] = 0
        }
    }
    
    #Case 2
    if (length(zeroToOne) == 0 & length(oneToZero) != 0){
        netSeq_path = list()
        #############
        # if (w2==1 & is.empty(y2)){
        #   w2 = 0
        # }
        #############
        if (w2 == 1){
            binSeq_path = rep(c(w1, 0, 0, 1), c(1, length(oneToZero), 1, 1))
            netSeq_path[[1]] = y1
            index = 1
            for (i in oneToZero){
                pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[i]
                netSeq_path[[index+1]] = netSeq_path[[index]]
                netSeq_path[[index+1]][i] = 0
                netSeq_path[[index+1]][pos1] = 0
                index = index+1
            }
            #randomly delete one edge then add
            pos = sample_modified(which(netSeq_path[[index]] == 1 & upper.tri(netSeq_path[[index]])), 1)
            pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[pos]
            netSeq_path[[index + 1]] = netSeq_path[[index]]
            netSeq_path[[index+1]][pos] = 0
            netSeq_path[[index+1]][pos1] = 0
            index = index+1
            netSeq_path[[index + 1]] = netSeq_path[[index-1]]
            
        }else{
            binSeq_path = rep(c(w1, 0), c(1, length(oneToZero)))
            netSeq_path[[1]] = y1
            index = 1
            for (i in oneToZero){
                pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[i]
                netSeq_path[[index+1]] = netSeq_path[[index]]
                netSeq_path[[index+1]][i] = 0
                netSeq_path[[index+1]][pos1] = 0
                index = index+1
            }
        }
    }
    
    #Case 3
    if(length(oneToZero) == 0 & length(zeroToOne) != 0){
        netSeq_path = list()
        #############
        # if (w2==0 & is.full(y2)){
        #   w2 = 1
        # }
        #############
        if (w2 == 0){
            binSeq_path = rep(c(w1, 1, 1, 0), c(1, length(zeroToOne), 1, 1))
            netSeq_path[[1]] = y1
            index = 1
            for (i in zeroToOne){
                pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[i]
                netSeq_path[[index+1]] = netSeq_path[[index]]
                netSeq_path[[index+1]][i] = 1
                netSeq_path[[index+1]][pos1] = 1
                index = index+1
            }
            #randomly add one edge then delete
            pos = sample_modified(which(netSeq_path[[index]] == 0 & upper.tri(netSeq_path[[index]])), 1)
            pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[pos]
            netSeq_path[[index + 1]] = netSeq_path[[index]]
            netSeq_path[[index+1]][pos] = 1
            netSeq_path[[index+1]][pos1] = 1
            index = index+1
            netSeq_path[[index + 1]] = netSeq_path[[index-1]]
            
        }else{
            binSeq_path = rep(c(w1, 1), c(1, length(zeroToOne)))
            netSeq_path[[1]] = y1
            index = 1
            for (i in zeroToOne){
                pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[i]
                netSeq_path[[index+1]] = netSeq_path[[index]]
                netSeq_path[[index+1]][i] = 1
                netSeq_path[[index+1]][pos1] = 1
                index = index+1
            }
        }
    }
    
    #Case 4
    if(length(oneToZero) != 0 & length(zeroToOne) != 0){
        netSeq_path = list()
        if (w2 == 0){
            binSeq_path = rep(c(w1, 1, 0), c(1, length(zeroToOne), length(oneToZero)))
            netSeq_path[[1]] = y1
            index = 1
            for (i in zeroToOne){
                pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[i]
                netSeq_path[[index+1]] = netSeq_path[[index]]
                netSeq_path[[index+1]][i] = 1
                netSeq_path[[index+1]][pos1] = 1
                index = index+1
            }
            for (i in oneToZero){
                pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[i]
                netSeq_path[[index+1]] = netSeq_path[[index]]
                netSeq_path[[index+1]][i] = 0
                netSeq_path[[index+1]][pos1] = 0
                index = index+1
            }
            
        }else{
            binSeq_path = rep(c(w1, 0, 1),c(1, length(oneToZero), length(zeroToOne)))
            netSeq_path[[1]] = y1
            index = 1
            for (i in oneToZero){
                pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[i]
                netSeq_path[[index+1]] = netSeq_path[[index]]
                netSeq_path[[index+1]][i] = 0
                netSeq_path[[index+1]][pos1] = 0
                index = index+1
            }
            for (i in zeroToOne){
                pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[i]
                netSeq_path[[index+1]] = netSeq_path[[index]]
                netSeq_path[[index+1]][i] = 1
                netSeq_path[[index+1]][pos1] = 1
                index = index+1
            }
        }
    }
    
    if (length(netSeq_path) == 2){
        netSeq_path[[4]] = netSeq_path[[2]]
        netSeq_path[[3]] = netSeq_path[[1]]
        netSeq_path[[2]] = netSeq_path[[1]]
        if (is.empty(netSeq_path[[1]])){
            nonEdge = which(netSeq_path[[1]] == 0 & upper.tri(netSeq_path[[1]]))
            pos = sample_modified(nonEdge,1)
            pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[pos]
            netSeq_path[[2]][pos] = 1
            netSeq_path[[2]][pos1] = 1
            binSeq_path = c(binSeq_path[1], 1, 0, binSeq_path[2])
        }else{
            edge = which(netSeq_path[[1]] == 1 & upper.tri(netSeq_path[[1]]))
            pos = sample_modified(edge, 1)
            pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[pos]
            netSeq_path[[2]][pos] = 0
            netSeq_path[[2]][pos1] = 0
            binSeq_path = c(binSeq_path[1], 0, 1, binSeq_path[2])
        }
        # binSeq_path = isPath_GetBinaryVar(netSeq_path)$binSeq_path
    }
    return (list(netSeq_path = netSeq_path, binSeq_path = binSeq_path))
}

#get changed edge indices at each transitioning time from netSeq_path
get_edgeID_change = function(netSeq_path, N, length){
    id = rep(0, length)
    for(i in 2:length){
        id[i] = which(netSeq_path[[i]]-netSeq_path[[i-1]] != 0 & upper.tri(netSeq_path[[i]]))
        # id[i] = get_triangular(which(netSeq_path[[i]]-netSeq_path[[i-1]] !=0), N)
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
        edge = which(netSeq_path_prop[[i-1]] == 1 & upper.tri(netSeq_path_prop[[i-1]]))
        # edge = get_triangular(which(netSeq_path_prop[[i-1]] == 1),N)
        netSeq_path_prop[[i]] = netSeq_path_prop[[i-1]]
        # pos = get_pos(edgeID_change_prop[i], N)
        pos = edgeID_change_prop[i]
        pos1 = matrix(1:(N^2), N, N, byrow = TRUE)[pos]
        if (is.element(edgeID_change_prop[i], edge)){
            netSeq_path_prop[[i]][pos] = 0
            netSeq_path_prop[[i]][pos1] = 0
            binSeq_path_prop[i] = 0
        }else{
            netSeq_path_prop[[i]][pos] = 1
            netSeq_path_prop[[i]][pos1] = 1
            binSeq_path_prop[i] = 1
        }
    }
    return(list(netSeq_path_prop = netSeq_path_prop, binSeq_path_prop = binSeq_path_prop))
}


#
# insert_pair = function()
#
# delete_pair = function()
#
# permute = function()


# Get proposal path from the current path: netSeq_path, binSeq_path
# by 3 types of small changes
# Return: netSeq_path_prop, binSeq_path_prop
get_proposalPath2 = function(netSeq_path, binSeq_path, N){
    length = length(netSeq_path)
    edgeID_change = get_edgeID_change(netSeq_path, N, length)
    
    action = sample(c(1, 2, 3), 1)
    #netSeq_path_prop = list()
    
    #paired insertion
    if (action == 1){
        cand = sample(which(upper.tri(netSeq_path[[1]])), 1)
        position = sort(sample_modified(2:(length-1), 2))
        r1 = position[1]
        r2 = position[2]
        edgeID_change_prop = update_edgeID_insert(cand, r1, r2, edgeID_change, length)
        res = getPathFrom_EdgeID(edgeID_change_prop, netSeq_path, N)
        netSeq_path_prop = res$netSeq_path_prop
        binSeq_path_prop = res$binSeq_path_prop
        #isPath_GetBinaryVar(netSeq_path_prop)
        #cat("Action 1 ")
        # if(!isPath_GetBinaryVar(netSeq_path_prop)$isPath){
        #   stop("Get wrong path by minor change!!!")
        # }
    }
    
    #paired deletion
    if (action == 2){
        edgeID_change_part = edgeID_change[-length(edgeID_change)]
        if (length(edgeID_change_part) == length(unique(edgeID_change_part)) | length(edgeID_change)<=4){
            #if (length(edgeID_change_part) == length(unique(edgeID_change_part))){
            # no pair to delete, permute instead
            permutedChange = sample_modified(edgeID_change[2:(length-1)], length-2)
            edgeID_change_prop = c(edgeID_change[1], permutedChange, edgeID_change[length])
            res = getPathFrom_EdgeID(edgeID_change_prop, netSeq_path, N)
            netSeq_path_prop = res$netSeq_path_prop
            binSeq_path_prop = res$binSeq_path_prop
            
            # # no pair to delete, insert instead
            # cand = sample(get_triangular(1:N^2, N),1)
            # position = sort(sample_modified(2:(length-1),2))
            # r1 = position[1]
            # r2 = position[2]
            # edgeID_change_prop = update_edgeID_insert(cand, r1, r2, edgeID_change, length)
            # res = getPathFrom_EdgeID(edgeID_change_prop, netSeq_path, N)
            # netSeq_path_prop = res$netSeq_path_prop
            # binSeq_path_prop = res$binSeq_path_prop
            #cat("Action 2-1 ")
            # if(!isPath_GetBinaryVar(netSeq_path_prop)$isPath){
            #   stop("Get wrong path by minor change!!!")
            # }
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
            # if(!isPath_GetBinaryVar(netSeq_path_prop)$isPath){
            #   stop("Get wrong path by minor change!!!")
            # }
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
        # if(!isPath_GetBinaryVar(netSeq_path_prop)$isPath){
        #   stop("Get wrong path by minor change!!!")
        # }
    }
    
    return (list(netSeq_path_prop = netSeq_path_prop, binSeq_path_prop = binSeq_path_prop))
}



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
    N = nrow(y1)
    path_shortest = get_shortestPath(y1, w1, y2, w2)
    netSeq_path_cur = path_shortest$netSeq_path
    binSeq_path_cur = path_shortest$binSeq_path
    prop_num = 1
    reject_num = 0
    path_num = 1
    netSeq_path_D = array(list(), D)
    binSeq_path_D = list()
    while (path_num <= D){
        temp = get_proposalPath2(netSeq_path_cur, binSeq_path_cur, N)
        prop_num = prop_num + 1
        netSeq_path_prop = temp$netSeq_path_prop
        binSeq_path_prop = temp$binSeq_path_prop
        
        # if(!isPath_GetBinaryVar(netSeq_path_prop)$isPath){
        #   stop("Get wrong proposed path!!!")
        # }
        # accept = min(1,exp(getPathProb(netSeq_path_prop, binSeq_path_prop, t1,t2, p, q, gamma))/
        #                exp(getPathProb(netSeq_path_cur, binSeq_path_cur, t1,t2, p, q, gamma)))
        accept = min(0, getPathProb(netSeq_path_prop, binSeq_path_prop, t1,t2, p, q, gamma) -
        getPathProb(netSeq_path_cur, binSeq_path_cur, t1,t2, p, q, gamma) )
        ###########
        # if(is.nan(accept)){
        #   accept = 0
        # }
        ##########
        if (sample(c(1,0), 1, prob = c(exp(accept), 1-exp(accept))) == 1 | prop_num > prop_max){
            # if (sample(c(1,0), 1, prob = c(exp(accept), 1-exp(accept))) == 1){
            netSeq_path_D[[path_num]] = netSeq_path_prop
            binSeq_path_D[[path_num]] = binSeq_path_prop
            netSeq_path_cur = netSeq_path_prop
            binSeq_path_cur = binSeq_path_prop
            prop_num = 0
            path_num = path_num + 1
            # if(!isPath_GetBinaryVar(netSeq_path_cur)$isPath){
            #   stop("Get wrong path!!!")
            # }
        }else{
            reject_num = reject_num + 1
        }
    }
    rate_accept = D/(D + reject_num)
    #cat("Got D Paths!\n")
    return (list(netSeq_path_D = netSeq_path_D, binSeq_path_D = binSeq_path_D,
    rate_accept = rate_accept))
}

getSamplePath_nonparallel = function(m, D, Y_true, W_true, t_obs, p_cur, q_cur,
gamma_cur, process = "ER", prop_max, burnIn_mcmc = 0){
    cat(m, "..")
    M = length(t_obs)
    sum_01 = 0
    sum_0 = 0
    sum_10 = 0
    sum_1 = 0
    sum_R = 0
    MHsamples = getSamplePath(D, y1 = Y_true[[m]], w1 = W_true[m], y2 = Y_true[[m+1]], w2 = W_true[m+1],
    t1 = t_obs[m], t2 = t_obs[m+1], p_cur, q_cur, gamma_cur,
    process = "ER", prop_max = prop_max)
    for (d in (burnIn_mcmc+1):D){
        counts = getCount(MHsamples$binSeq_path_D[[d]])
        sum_01 = sum_01 + counts[1]
        sum_0 = sum_0 + counts[2]
        sum_10 = sum_10 + counts[3]
        sum_1 = sum_1 + counts[4]
        sum_R = sum_R + counts[5]
    }
    return (c(sum_01, sum_0, sum_10, sum_1, sum_R))
}

getSamplePath_parallel = function(processID, D, partic, t_obs, p_cur, q_cur,
gamma_cur, process = "ER", prop_max, ids, burnIn = 0, burnIn_mcmc = 0){
    source("MHSamples.R")
    M = length(t_obs)
    sum_01 = 0
    sum_0 = 0
    sum_10 = 0
    sum_1 = 0
    sum_R = 0
    m = processID %% (M-1-burnIn)
    h = processID %/% (M-1-burnIn) + 1
    if (m == 0){
        m = M - 1 - burnIn
        h = h - 1
        cat("  ", h, "...")
    }
    b = ids[h]
    Y_true = partic$particles[[b]]
    W_true = partic$W[b,]
    m = m + burnIn
    MHsamples = getSamplePath(D, y1 = Y_true[[m]], w1 = W_true[m], y2 = Y_true[[m+1]], w2 = W_true[m+1],
    t1 = t_obs[m], t2 = t_obs[m+1], p_cur, q_cur, gamma_cur,
    process = "ER", prop_max = prop_max)
    for (d in (burnIn_mcmc+1):D){
        counts = getCount(MHsamples$binSeq_path_D[[d]])
        sum_01 = sum_01 + counts[1]
        sum_0 = sum_0 + counts[2]
        sum_10 = sum_10 + counts[3]
        sum_1 = sum_1 + counts[4]
        sum_R = sum_R + counts[5]
    }
    return (c(sum_01, sum_0, sum_10, sum_1, sum_R))
}

# getSamplePath_parallel = function(processID, D, Y_true, W_true, t_obs, p_cur, q_cur,
#                                   gamma_cur, process = "ER", prop_max){
#   source("MHSamples.R")
#   M = length(t_obs)
#   sum_01 = 0
#   sum_0 = 0
#   sum_10 = 0
#   sum_1 = 0
#   sum_R = 0
#   m = processID + 1
#   MHsamples = getSamplePath(D, y1 = Y_true[[m-1]], w1 = W_true[m-1], y2 = Y_true[[m]], w2 = W_true[m],
#                             t1 = t_obs[m-1], t2 = t_obs[m], p_cur, q_cur, gamma_cur,
#                             process = "ER", prop_max = prop_max)
#   for (d in 1:D){
#     counts = getCount(MHsamples$binSeq_path_D[[d]])
#     sum_01 = sum_01 + counts[1]
#     sum_0 = sum_0 + counts[2]
#     sum_10 = sum_10 + counts[3]
#     sum_1 = sum_1 + counts[4]
#     sum_R = sum_R + counts[5]
#   }
#   return (c(sum_01, sum_0, sum_10, sum_1, sum_R))
# }


# get count statistics from a series of MH samples
getCount = function(binSeq_path_D){
    len = length(binSeq_path_D)
    # sum_R = len - 2
    sum_R = len - 1
    sum_1 = sum(binSeq_path_D[1:(len-1)])
    sum_0 = len - 1 - sum_1
    diff = binSeq_path_D[1:(len-1)] - binSeq_path_D[2:len]
    sum_10 = sum(diff == 1)
    sum_01 = sum(diff == (-1))
    return (c(sum_01, sum_0, sum_10, sum_1, sum_R))
}



