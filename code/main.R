source("GetSimSeq.R")
source("EM-parallel.R")

set.seed(1022)
M = 10
N = 10
p = 0.7; q = 0.3; gamma = 2; alpha = 0.3; beta = 0.2

#random draw M time points uniformly from (0, 20)
t_obs = sort(runif(M, 0, 20))
#or
#t_obs = c(0.5, 1.3, 4, 4.5, 8, 10)

#initialize the starting network with a random graph with 1/5 of edges
net_start = matrix(0, N, N)
edges = sample_modified(get_triangular(1:N^2, N), ceiling(choose(N,2)/5))
for (i in 1:length(edges)){
    pos = get_pos(edges[i], N)
    net_start[pos[1], pos[2]] = 1
    net_start[pos[2], pos[1]] = 1
}

res = get_SimSeq (p, q, gamma, N, t_obs, net_start = net_start)

net_obs = res$network_seq_obs
net_obs_noise = add_Noise(net_obs, alpha, beta)
t_all = res$t_transition
net_hidden_all = res$network_seq_hidden

#t_obs # observation moments
#net_obs_noise #observed networks with noise

#net_obs # noise-free networks at observation moments
#t_all #all the transitioning time points in the hidden layer
#net_hidden_all #all the true networks at transitioning time points in the hidden layer

set.seed(1022)
Y = net_obs_noise

##############################################
init = c(0.7, 0.3, 2, 0.3, 0.2)
B = 50000
thr = 1e-1
D = 200
H = 1000
burnIn = ceiling(0.2*M); phi = B*0.9
MaxIter = 30
prop_max = 40
particleParallel = F; clusterSizePF = 3
MHParallel = F; clusterSizeMH = 3
outFile = "output.txt"


cat("Settings: p, q, gamma, alpha, beta = ", c(p, q, gamma, alpha, beta), "\n")
cat("Settings: p, q, gamma, alpha, beta = ", c(p, q, gamma, alpha, beta), "\n", sep = " ", file = outFile, append = TRUE)
cat("          N, M, B, H, D", c(N, M, B, H, D), "\n")
cat("          N, M, B, H, D", c(N, M, B, H, D), "\n", sep = " ", file = outFile, append = TRUE)

time = proc.time()
results = EM_parallel (Y, t_obs, B, process = "ER", thr, init = init, D = D, MaxIter = MaxIter, 
                       prop_max = prop_max, burnIn, phi, H, particleParallel, MHParallel, outFile = "output.txt")
time_used = proc.time() - time
time_used

cat("Time elapsed:", time_used[3],"\n")
cat("@@@@@@@RESULTS@@@@@@@:\n")
cat("Initialization:", init, "\n")
cat("Truth:", c(p, q, gamma, alpha, beta), "\n")
cat("Settings: (N, M, B, H, D)", c(N, M, B, H, D), "\n")
cat("(thr, burnIn, phi, MaxIter):", c(thr, burnIn, phi, MaxIter), "\n")
cat("(particleParallel, clusterSizePF):", c(particleParallel, clusterSizePF), "\n")
cat("Estimation(p,q,gamma,alpha,beta): \n", results$Estimation, "\n")
cat("iteration:\n", results$iteration, "\n")


cat("Time elapsed:", time_used[3],"\n", sep = " ", file = outFile, append = TRUE)
cat("@@@@@@@RESULTS@@@@@@@:\n", sep = " ", file = outFile, append = TRUE)
cat("Initialization:", init, "\n", sep = " ", file = outFile, append = TRUE)
cat("Truth:", c(p, q, gamma, alpha, beta), "\n", sep = " ", file = outFile, append = TRUE)
cat("Settings: (N, M, B, H, D)", c(N, M, B, H, D), "\n", sep = " ", file = outFile, append = TRUE)
cat("(thr, burnIn, phi, MaxIter):", c(thr, burnIn, phi, MaxIter), "\n", sep = " ", file = outFile, append = TRUE)
cat("(particleParallel, clusterSizePF):", c(particleParallel, clusterSizePF), "\n", sep = " ", file = outFile, append = TRUE)
cat("Estimation(p,q,gamma,alpha,beta): \n", results$Estimation, "\n", sep = " ", file = outFile, append = TRUE)
cat("iteration:\n", results$iteration, "\n", sep = " ", file = outFile, append = TRUE)
