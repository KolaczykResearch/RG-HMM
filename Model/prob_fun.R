source("utilities.R")
##########################################
###This is for PARTICLE FILTERING#########
##########################################
# get log sum of two exponentials, ie get l from l1, l2 where 
# l = log (exp(l1) + exp(l2)) = l1 + log(1 + exp(l2 - l1)) if l1>l2 
# l1 = log(p1), l2 = log(p2), l = log(p1+p2)
lse2 = function (l1, l2){
  if (l1==-Inf & l2==-Inf){
    return(-Inf)
  }else{
    m = max(l1, l2)
    d = -abs(l1-l2)
    return (ifelse(d < log(.Machine$double.eps), m, m + log(1 + exp(d))))
  }
}

# get log sum of k exponentials, ie get l from l1, l2,..., lk
# l = log(exp(l1) + ...+ exp(lk))
# input: ls, which is vector of [l1, l2, ..., lk]
lse = function(ls){
  Reduce (lse2, ls)
}

#the conditional prob. function of observed network given the true (log)
#input: two matrices, type I/II errors
#output: numeric value
#getCondProb = function (observe, true, alpha, beta){
#    countVec = getCounts_abcd (observe, true)
#    a = countVec[1]
#    b = countVec[2]
#    c = countVec[3]
#    d = countVec[4]
#    loglik = c*log(alpha) + d*log(1-alpha)+ b*log(beta) + a*log(1-beta)
#    return (loglik)
#}
getCondProb = function (observe, true, alpha, beta){
    countVec = getCounts_abcd (observe, true)
    a = countVec[1]
    b = countVec[2]
    c = countVec[3]
    d = countVec[4]
    if(alpha == 0 | beta == 0){
      loglik = log (alpha^c * (1-alpha)^d * beta^b * (1-beta)^a)
    }else{
      loglik = c*log(alpha) + d*log(1-alpha)+ b*log(beta) + a*log(1-beta)
    }
    return (loglik)
}
#get a, b, c, d (counts corresponding to type I and type II errors)
#based on observe and true
#input: two matrices
#output: array of 4 integers
getCounts_abcd = function (observe, true){
    N = nrow(observe)
    diffMat = observe - true
    a = length(which(true == 1 & observe == 1))/2
    b = length(which((diffMat)==-1))/2 #1->0
    c = length(which((diffMat)==1))/2 #0->1
    d = (length(which(true == 0 & observe == 0)) - N)/2
    return (c(a,b,c,d))
}

#observe = net_obs_noise[[3]]
#true = net_obs[[3]]

##########################################
###This is for METROPOLIS HASTINGS SAMPLES#########
##########################################
#transition probability function log_f(y2,w2|y1,w1) = log_g(y2|w2,y1) + log_h(w2|y1,w1)

#Conditional log pdf of y2 given w2, y1
# log_g = function(y2,w2,y1){
#   res = checkAdding(y1,y2)
#   if (res == 1 & w2==1){
#     lik = 1/nonedge_num(y1)
#   }else if(res == 0 & w2==0){
#     lik = 1/edge_num(y1)
#   }else{
#     lik = 0
#   }
#   return (log(lik))
# }
# assume feasible, log pdf of y2 given w2, y1 for ER process
log_g_ER = function(y2, w2, y1){
    if (w2 == 1){
        loglik = - log(nonedge_num(y1))
    }else{
        loglik = -log(edge_num(y1))
    }
    return (loglik)
}
# assume feasible, log pdf of y2 given w2, y1 for PR process
log_g_PR = function(y2, w2, y1){
  graph = graph_from_adjacency_matrix(y1, mode = "undirected")
  clusters = clusters(graph)
  if (clusters$no == 1){
    loglik = log_g_ER(y2, w2, y1)
  }else if (w2 == 1){
    # add edge
    ei = which((y2 - y1) == 1 & upper.tri(y1))
    prod_score_add = get_prod_score_add(y1, clusters)
    nonedge = nrow(prod_score_add)
    if (nonedge == 1){
      loglik = log(1)
    }else{
      score_i = prod_score_add[prod_score_add[, 1]==ei, 2]
      non_i_bool = (prod_score_add[, 1]!=ei)
      greater = sum(non_i_bool & prod_score_add[, 2]>score_i)
      equal = sum(non_i_bool & prod_score_add[, 2]==score_i)
      loglik = log(greater + equal/2) - log(choose(nonedge, 2))
    }
  }else{
    # delete edge
    ei = which((y2 - y1) == -1 & upper.tri(y1))
    prod_score_delete = get_prod_score_delete(y1, clusters)
    edge = nrow(prod_score_delete)
    if (edge == 1){
      loglik = log(1)
    }else{
      score_i = prod_score_delete[prod_score_delete[, 1]==ei, 2]
      non_i_bool = (prod_score_delete[, 1]!=ei)
      less = sum(non_i_bool & prod_score_delete[, 2]<score_i)
      equal = sum(non_i_bool & prod_score_delete[, 2]==score_i)
      loglik = log(less + equal/2) - log(choose(edge, 2))
    }
  }
  return (loglik)
}

#Conditional log pdf of w2 given w1, y1
log_h = function(w2,y1,w1,p,q){
    if (is.empty(y1)){
        loglik = log(w2)
    }else if (is.full(y1)){
        loglik = log(1-w2)
    }else{
        if (w1 == 0){
            # loglik = w2*log(p)+(1-w2)*log(1-p)
            loglik = log(p^w2 * (1-p)^(1-w2))
        }
        if (w1 == 1){
            # loglik = w2*log(1-q)+(1-w2)*log(q)
            loglik = log((1-q)^w2 * q^(1-w2))
        }
    }
    return (loglik)
}

#y1 = netSeq_path[[4]]
#y2 = netSeq_path[[5]]
#w1 = binSeq_path[[4]]
#w2 = binSeq_path[[5]]

#the probability function (log) for feasible sample path (S,V) that brings y(t1) to y(t2),
#conditional on Y(t1), W(t2)
#input: netSeq_path from y1 to y2, binSeq_path from w1 to w2, parameters
#t1, t2: obs. moments for y1 and y2, respectively
#output: loglik
#remark: s1 = netSeq_path[-1], v1 = binSeq_path[-1], y1 = netSeq_path[1], y2 = netSeq_path[lastItem] (the last change is y2)
#example path that brings y(t1) to y(t2) is: [y(t1), s1, s2, s3, y(t2)] (3 changes); Since last change s3=y(t2), then netSeq_path is of the form [y(t1), s1, s2, y(t3)]
#if length of netSeq_path is 1, it means that y(t1)=y(t2), there is no change between t1 and t2
getPathProb =  function (netSeq_path, binSeq_path, t1,t2, p, q, gamma, process = "ER"){
    M = length(netSeq_path)
    s1 = netSeq_path[-1] # exclude the first one
    v1 = binSeq_path[-1]
    R1 = M-1
    # loglik = (-1)*gamma*(t2-t1)+R1*log(gamma*(t2-t1))-log(factorial(R1))
    loglik = dpois(R1, gamma*(t2-t1), log = T)
    if (R1 == 0){
      return(loglik) # there is no change between t1, t2
    }
    if(process == "ER"){
      for (r in 2:(R1+1)){
        loglik = loglik + log_g_ER(y2=netSeq_path[[r]], w2=binSeq_path[[r]], y1=netSeq_path[[r-1]]) +
          log_h(w2=binSeq_path[[r]],y1=netSeq_path[[r-1]],w1=binSeq_path[[r-1]],p,q)
      }
    }
    if(process == "PR"){
      for (r in 2:(R1+1)){
        loglik = loglik + log_g_PR(y2=netSeq_path[[r]], w2=binSeq_path[[r]], y1=netSeq_path[[r-1]]) +
          log_h(w2=binSeq_path[[r]],y1=netSeq_path[[r-1]],w1=binSeq_path[[r-1]],p,q)
      }
    }
    #######
    # if (is.nan(loglik)){
    #   loglik = -Inf
    # }
    #######
    return(loglik)
}

#netSeq_path = net_hidden_all[1:7]
#binSeq_path = isPath_GetBinaryVar(netSeq_path)$binSeq_path
#binSeq_path
#t1 = t_all[1]
#t2 = t_all[3]
#getPathProb(netSeq_path, binSeq_path, t1,t2, p, q, gamma)

#derivative of h wrt p,q
h_prime_p = function(w2,y1,w1){
    if (is.empty(y1)){
        deriv = 0
    }else if (is.full(y1)){
        deriv = 0
    }else{
        if (w2==0){
            deriv = -1
        }
        if (w2 ==1){
            deriv = 1
        }
    }
    return (deriv)
}

h_prime_q = function(w2,y1,w1){
    if (is.empty(y1)){
        deriv = 0
    }else if (is.full(y1)){
        deriv = 0
    }else{
        if (w2==0){
            deriv = 1
        }
        if (w2 ==1){
            deriv = -1
        }
    }
    return (deriv)
}



#get score_m
#input: netSeq_path from y1 to y2, binSeq_path from w1 to w2, parameters
#t1, t2: obs. moments for y1 and y2, respectively
#output: score function wrt p,q,gamma
#remark: s1 = netSeq_path[-1], v1 = binSeq_path[-1], y1 = netSeq_path[1], y2 = netSeq_path[lastItem]
get_score_m = function(netSeq_path, binSeq_path, t1,t2, p, q, gamma){
    M = length(netSeq_path)
    R1 = M-1
    score_m_p = 0
    score_m_q = 0
    score_m_gamma = -(t2-t1) + R1/gamma
    
    for (r in 2:(R1+1)){
        score_m_p = score_m_p + 1/exp(log_h(w2=binSeq_path[[r]],
        y1=netSeq_path[[r-1]], w1=binSeq_path[[r-1]],p,q))*
        h_prime_p(w2=binSeq_path[[r]],
        y1=netSeq_path[[r-1]], w1=binSeq_path[[r-1]])
        score_m_q = score_m_q + 1/exp(log_h(w2=binSeq_path[[r]],
        y1=netSeq_path[[r-1]], w1=binSeq_path[[r-1]],p,q))*
        h_prime_q(w2=binSeq_path[[r]],
        y1=netSeq_path[[r-1]], w1=binSeq_path[[r-1]])
        #cat(score_m_p,score_m_q, "\n")
    }
    return (c(score_m_p, score_m_q, score_m_gamma))
}

#get gradient of score_m
#input: netSeq_path from y1 to y2, binSeq_path from w1 to w2, parameters
#t1, t2: obs. moments for y1 and y2, respectively
#output: derivative of score function wrt p, q, gamma
#remark: s1 = netSeq_path[-1], v1 = binSeq_path[-1], y1 = netSeq_path[1], y2 = netSeq_path[lastItem]
get_score_m_grad = function(netSeq_path, binSeq_path, t1,t2, p, q, gamma){
    M = length(netSeq_path)
    R1 = M-1
    score_m_grad_p = 0
    score_m_grad_q = 0
    score_m_grad_gamma = - R1/(gamma)^2
    
    for (r in 2:(R1+1)){
        score_m_grad_p = score_m_grad_q - 1/(exp(log_h(w2=binSeq_path[[r]],
        y1=netSeq_path[[r-1]], w1=binSeq_path[[r-1]],p,q)))^2*
        (h_prime_p(w2=binSeq_path[[r]], y1=netSeq_path[[r-1]], w1=binSeq_path[[r-1]]))^2
        
        score_m_grad_q = score_m_grad_q - 1/(exp(log_h(w2=binSeq_path[[r]],
        y1=netSeq_path[[r-1]], w1=binSeq_path[[r-1]],p,q)))^2*
        (h_prime_q(w2=binSeq_path[[r]], y1=netSeq_path[[r-1]], w1=binSeq_path[[r-1]]))^2
        #cat(score_m_p,score_m_q, "\n")
    }
    return (c(score_m_grad_p, score_m_grad_q, score_m_grad_gamma))
}


#get score_m and score_m_grad wrt the natural parameters
#input: netSeq_path from y1 to y2, binSeq_path from w1 to w2, parameters
#t1, t2: obs. moments for y1 and y2, respectively
#output: score function wrt p_til = p_til(p), q_til = q_til(q),gamma_til = log(gamma)
#remark: s1 = netSeq_path[-1], v1 = binSeq_path[-1], y1 = netSeq_path[1], y2 = netSeq_path[lastItem]
get_score_m_natural = function(netSeq_path, binSeq_path, t1,t2, p, q, gamma){
    score_m_natural = rep(0,3)
    score_m_grad_natural = rep(0,3)
    p_til = log(p/(1-p))
    q_til = log(q/(1-q))
    gamma_til = log(gamma)
    score_m = get_score_m(netSeq_path, binSeq_path, t1,t2, p, q, gamma)
    score_m_grad = get_score_m_grad(netSeq_path, binSeq_path, t1,t2, p, q, gamma)
    score_m_natural[1] = score_m[1]*sigmoid_dev(p_til)
    score_m_natural[2] = score_m[2]*sigmoid_dev(q_til)
    score_m_natural[3] = score_m[3]*exp(gamma_til)
    score_m_grad_natural[1] = score_m_grad[1]*(sigmoid_dev(p_til))^2 + score_m[1]*sigmoid_dev_sec(p_til)
    score_m_grad_natural[2] = score_m_grad[2]*(sigmoid_dev(q_til))^2 + score_m[2]*sigmoid_dev_sec(q_til)
    score_m_grad_natural[3] = score_m_grad[3]*(exp(gamma_til))^2 + score_m[3]*exp(gamma_til)
    return (c(score_m_natural, score_m_grad_natural))
}



##get expected observed data score function
#getExpObdScore_obs = function()


##get expected log-likelihood

