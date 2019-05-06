source("DataSet.R")
source("EMParallel2.R")

# ===============================================
# Set parameters
# ===============================================
set.seed(10221211)
M = 10
N = 20
p = 0.7; q = 0.3; gamma = 2; alpha = 0.03; beta = 0.02
rate_obs = 0.6

# ===============================================
# Generate synthetic data
# ===============================================
param_val = c(p, q, gamma, alpha, beta)
data = DataSet (N, M, param_val, rate_obs, process = "ER")

# Data observed
t_obs = data$t_obs # observation moments
net_obs_noise = data$net_obs_noise # observed networks with noise

# Ground truth in hidden layer
net_obs = data$net_obs # true networks at observational moments 
bin_obs = data$bin_obs # binary variables at observational moments
net_hidden = data$net_hidden # true networks at all transitioning moments 
bin_hidden = data$bin_hidden # binary variables at all transitioning moments
t_all = data$t_all # all transitioning time points 

# ===============================================
# Run EM on synthetic data
# ===============================================
set.seed(10221211)
Y = net_obs_noise
init = c(p, q, gamma, alpha, beta)
#init = c(0.1, 0.1, 1, 0.01, 0.01)
B = 5000
thr = 1e-1
D = 10
H = 10
burnIn = ceiling(0.0*M); phi = B*0.9
burnIn_mcmc = ceiling(0.0*D)
MaxIter = 30
prop_max = 50
particleParallel = F; # clusterSizePF = 3
MHParallel = F # clusterSizeMH = 3
process = 'ER'
#######################
# generating particles
#######################
p_cur = init[1]
q_cur = init[2]
gamma_cur = init[3]
alpha_cur = init[4]
beta_cur = init[5]
cat("Initialization:", init, "\n")
M = length(t_obs)
# stopping condition for EM:
stopConEM = 1
iter = 0
# Particle Filtering
time1 = proc.time()
partic = particleFil(B, Y, t_obs, p_cur, q_cur, gamma_cur, alpha_cur, beta_cur,
                      process = process, E_step = F)
time_used2 = proc.time() - time1
cat("  Time elapsed:", time_used2[3],"\n")

#######################
# update alpha beta
#######################
B = length(partic$particles)
M = length(partic$particles[[1]])
N = nrow(partic$particles[[1]][[1]])
# set = sample(B, phi)
# sample phi last indices
set = sample(1:B, phi, replace = F, prob = partic$weights)
#set = which.max(partic$weights) 
sum_c = 0
sum_cd = 0
sum_b = 0
sum_ab = 0
# ==
c = matrix(0, phi, M)
cd = matrix(0, phi, M)
b = matrix(0, phi, M)
ab = matrix(0, phi, M)
id = 0
# ==
for (i in set){
  id = id + 1
  temp = get_particle_path(partic, i, M)
  Y_true = temp$Y_true
  for (j in (burnIn+1):M){
    sum_c = sum_c + sum((Y[[j]] - Y_true[[j]])==1)/2
    sum_cd = sum_cd + nonedge_num(Y_true[[j]])
    sum_b = sum_b + sum((Y[[j]] - Y_true[[j]])==-1)/2
    sum_ab = sum_ab + edge_num(Y_true[[j]])
    # ==
    c[id, j] = sum((Y[[j]] - Y_true[[j]])==1)/2
    cd[id, j] = nonedge_num(Y_true[[j]])
    b[id, j] = sum((Y[[j]] - Y_true[[j]])==-1)/2
    ab[id, j] = edge_num(Y_true[[j]])
    # ==
  }
}
cat(c(sum_c, sum_cd, sum_b, sum_ab)/(phi*M))
alpha_hat = sum_c/sum_cd
beta_hat = sum_b/sum_ab
cat("alpha = ", alpha_hat, ", beta = ", beta_hat)
  
hist(apply(c, 1, sum))
hist(apply(cd, 1, sum))
hist(apply(b, 1, sum))
hist(apply(ab, 1, sum))
hist(apply(c, 2, sum)/phi)
hist(apply(cd, 2, sum)/phi)
hist(apply(b, 2, sum)/phi)
hist(apply(ab, 2, sum)/phi)

apply(c, 2, sum)/phi
apply(cd, 2, sum)/phi
apply(b, 2, sum)/phi
apply(ab, 2, sum)/phi

apply(c, 2, sum)/apply(cd, 2, sum)
apply(b, 2, sum)/apply(ab, 2, sum)
################################
# check alpha, beta with Y_true
##############################
Y_true = net_obs
Y = net_obs_noise
# check alpha, beta with Y_true
sum_c_t = 0
sum_cd_t = 0
sum_b_t = 0
sum_ab_t = 0
for (j in 1:M){
  sum_c_t = sum_c_t + sum((Y[[j]] - Y_true[[j]])==1)/2
  sum_cd_t = sum_cd_t + nonedge_num(Y_true[[j]])
  sum_b_t = sum_b_t + sum((Y[[j]] - Y_true[[j]])==-1)/2
  sum_ab_t = sum_ab_t + edge_num(Y_true[[j]])
}
alpha_hat = sum_c_t/sum_cd_t
beta_hat = sum_b_t/sum_ab_t
cat("alpha = ", alpha_hat, ", beta = ", beta_hat)
c(sum_c_t, sum_cd_t, sum_b_t, sum_ab_t)/M


################################
## Check distribution
#############################
# check alpha, beta with Y_true
Y_true = net_obs
Y = net_obs_noise
# check alpha, beta with Y_true
c_t = numeric(M)
cd_t = numeric(M)
b_t = numeric(M)
ab_t = numeric(M)
for (j in 1:M){
  c_t[j] = sum((Y[[j]] - Y_true[[j]])==1)/2
  cd_t[j] = nonedge_num(Y_true[[j]])
  b_t[j] = sum((Y[[j]] - Y_true[[j]])==-1)/2
  ab_t[j] = edge_num(Y_true[[j]])
}
alpha_hat = sum(c_t)/sum(cd_t)
beta_hat = sum(b_t)/sum(ab_t)
cat("alpha = ", alpha_hat, ", beta = ", beta_hat)
hist(c_t)
hist(cd_t)
hist(b_t) 
hist(ab_t)
c_t/cd_t
b_t/ab_t
