library(R.matlab) 
library(igraph)
source("auto_segment_finder_new.R")

# ---------------------------------------------
# Define helper functions to load and plot data
# ---------------------------------------------
# Get network time series Y and time sequence t_obs from MatLab data, where
## t_obs = data$nets[,,1]$t
## Y[[i]] = data$nets[,,1]$C[,,i]
get_Y_t = function(data){
  t_obs = data$nets[,,1]$t
  Y = list()
  for (i in 1:length(t_obs)){
    Y[[i]] = data$nets[,,1]$C[,,i]
  }
  return(list(Y=Y, t_obs=t_obs))
}

makeplot = function(data, range, title, metric='gcc'){
  # get gcc over time and t_all
  N = dim(data$nets[,,1]$C[,,1])[1]
  t_all = data$nets[,,1]$t
  gcc = numeric(0)
  if(metric=='gcc'){
    for (t in 1:length(t_all)){
      net_t = data$nets[,,1]$C[,,t]
      gcc = c(gcc, getSizeGCC(net_t, order=1)/N)
    }
  }else{
    for (t in 1:length(t_all)){
      net_t = data$nets[,,1]$C[,,t]
      gcc = c(gcc, get_density(net_t, N))
    }
  }
  if(metric=='gcc'){
    plot(t_all, gcc, type='l', 
         xlab='Time [s]', ylab='Proportion of nodes in GCC', 
         main = title, ylim=c(0,1))
  }else{
    plot(t_all, gcc, type='l', 
         xlab='Time [s]', ylab='Density', 
         main = title)
  }
  
  start_ = range[1]
  end_ = range[2]
  abline(v=start_, col='red', lwd=1, lty=2)
  abline(v=end_, col='blue', lwd=1, lty=2)
  return(list(gcc=gcc, t_all=t_all))
}

add_window = function(xleft, xright, col=rgb(1, 1, 0, 0.2)){
  rect(xleft, xright, 
       ybottom=-4, ytop=4, col=col)
}

color_segments = function(temp_plot, segments, level_low, level_high){
  cbbPalette = rep(c("#E69F00", "#56B4E9", "#D55E00", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#00AFBB"), 3)
  mat = segments
  for (i in 1:nrow(mat)){
    left = mat[i, 1]
    right = mat[i, 2]
    lines(temp_plot$t_all[left:right], temp_plot$gcc[left:right], col=cbbPalette[i], lwd=2)
  }
  abline(h=level_low, col='grey')
  abline(h=level_high, col='grey')
}

# ---------------------------------------------
# Load one set of data: S1, S2 or S3; and one ROI
# ---------------------------------------------
# S1 ===============================
# Left
data1_1 = readMat('../../DataSeizure/networks_fdr/nets_EP001_clip_1_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat')
t_start = 64
title = 'P1S1_1_0.5_left'
range = c(64, 314) # from seizure onset to termination 

# ROI is obtained from dynamic community tracking, see main_track_communities.m
ROI = c(62, 144)  # [1,]  126  134 
ROI = c(226, 250) # [1,]  453  463
ROI = c(304, 315) # [1,]  609  624


# S2 ===============================
# Left
data1_1 = readMat('../../DataSeizure/networks_fdr/nets_EP001_clip_2_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat')
t_start = 203
title = 'P1S2_1_0.5_left'
range = c(203, 412)
ROI = c(200, 230) # [1,]  236  253, 
ROI = c(230, 307) # [1,]  295  318
ROI = c(396, 410) # [1,]  627  646 

# S3 ===============================
# S3 left
data1_1 = readMat('../../DataSeizure/networks_fdr/nets_EP001_clip_3_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat')
t_start = 80
title = 'P1S3_1_0.5_left'
range = c(80, 250)
ROI = c(114, 126) # [1,]  229  238
ROI = c(126, 150) # [1,]  256  269
ROI = c(154, 190) # [1,]  309  316
ROI = c(191, 230) # [1,]  385  395

# ---------------------------------------------
# Get network time series: Y
#     observational time points: t_obs
# ---------------------------------------------
# get t_obs, Y from matlab data
temp = get_Y_t(data1_1)
t_obs = temp$t_obs
Y = temp$Y

# ---------------------------------------------
# Run auto_segment_finder_fine_tune()
# ---------------------------------------------
res = auto_segment_finder_fine_tune(Y, t_obs, ROI=ROI)
res$segments_fine_tuned
res$segments_all

# ---------------------------------------------
# Plot ramp-up(s)
# ---------------------------------------------
# pdf: 13X5.91
par(mfrow=c(2,1), mai = c(1, 1, 0.5, 0.1))

# Plots: after fine-tune
# gcc plots
temp_plot = makeplot(data1_1, range, title, metric='gcc') 
add_window(ROI[1], ROI[2])
color_segments(temp_plot, res$segments_fine_tuned, res$level_low_gcc, res$level_high_gcc)
# density plots
temp_plot = makeplot(data1_1, range, title, metric='den')
add_window(ROI[1], ROI[2])
color_segments(temp_plot, res$segments_fine_tuned, res$level_low_den, res$level_high_den)

# Plots: before fine-tune
# gcc plots
temp_plot = makeplot(data1_1, range, title, metric='gcc') 
add_window(ROI[1], ROI[2])
color_segments(temp_plot, res$segments_all, res$level_low_gcc, res$level_high_gcc)
# density plots
temp_plot = makeplot(data1_1, range, title, metric='den')
add_window(ROI[1], ROI[2])
color_segments(temp_plot, res$segments_all, res$level_low_den, res$level_high_den)



