library(R.matlab)
library('igraph')
source("auto_segment_finder3.R")

# ------------------------------
# Define plot related functions
# ------------------------------
# matlab data: 
# t_obs = data$nets[,,1]$t
# Y[[i]] = data$nets[,,1]$C[,,i]
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
         main = title, ylim=c(0,1))
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

color_segments_all = function(temp_plot, temp_segment){
  cbbPalette = rep(c("#E69F00", "#56B4E9", "#D55E00", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#00AFBB"), 3)
  splits = temp_segment$segments_all
  for (i in 2:length(splits)){
    lines(temp_plot$t_all[splits[i-1]:splits[i]], temp_plot$gcc[splits[i-1]:splits[i]], col=cbbPalette[i-1], lwd=2)
  }
  abline(h=temp_segment$level_low, col='grey')
  abline(h=temp_segment$level_high, col='grey')
}

color_segments_mat = function(temp_plot, temp_segment, level_low, level_high){
  cbbPalette = rep(c("#E69F00", "#56B4E9", "#D55E00", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#00AFBB"), 3)
  mat = temp_segment$segments_mat
  for (i in 1:nrow(mat)){
    left = mat[i, 1]
    right = mat[i, 2]
    lines(temp_plot$t_all[left:right], temp_plot$gcc[left:right], col=cbbPalette[i], lwd=2)
  }
  abline(h=level_low, col='grey')
  abline(h=level_high, col='grey')
}

# ------------------------------
# Load network data 
# ------------------------------
# S1 ===============================
# Left
data1_1 = readMat('Data/nets_EP001_clip_1_ecog_ws_1_wstp_0.5_sub_hemi_left.mat')
t_start = 64
title = 'P1S1_1_0.5_sub_Left'
range = c(64, 314)
# t_start = 250
# range = c(250, 314)

# S2 ===============================
# Left
data1_1 = readMat('Data/nets_EP001_clip_2_ecog_ws_1_wstp_0.5_sub_hemi_left.mat')
t_start = 203
title = 'P1S2_1_0.5_sub_Left'
range = c(203, 412)
# range = c(355, 412)

# S3 ===============================
# Left
data1_1 = readMat('Data/nets_EP001_clip_3_ecog_ws_1_wstp_0.5_sub_hemi_left.mat')
t_start = 80
title = 'P1S3_1_0.5_sub_Left'
range = c(80, 250)


# get t_obs, Y from matlab data
temp = get_Y_t(data1_1)
t_obs = temp$t_obs
Y = temp$Y

par(mfrow=c(2,1), mai = c(1, 1, 0.5, 0.1))
# pdf: 13X5.91

# ------------------------------
# Detect the proper curve to test with auto_segment_finder()
# ------------------------------
# 1. Firstly, detect segments based on gcc ------------------------------
temp_segment = auto_segment_finder(Y, t_obs, t_start, metric='gcc')
temp_plot = makeplot(data1_1, range, title, metric='gcc') # before tail off
level_low_gcc = temp_segment$level_low
level_high_gcc = temp_segment$level_high
gcc = temp_plot$gcc

# 2. Secondly, detect segments based on gcc ------------------------------
temp_segment = auto_segment_finder(Y, t_obs, t_start, metric='density')
temp_plot = makeplot(data1_1, range, title, metric='density') # before tail off
level_low_den = temp_segment$level_low
level_high_den = temp_segment$level_high
den = temp_plot$gcc

# 3. Lastly, take overlapping --------------------------------------------
temp_segment_overlapped = find_qualified_segments(Y, t_obs, t_start)
temp_segment_final = final_segments(temp_segment_overlapped$segments_mat, gcc, den, 
                                    level_low_gcc, level_high_gcc, level_low_den, level_high_den, t_obs, t_start)
temp_segment_final$segments_mat


# ------------------------------
# Make plots
# ------------------------------
#gcc
temp_plot = makeplot(data1_1, range, title, metric='gcc')
add_window(t_start-50, t_start+50)
color_segments_mat(temp_plot, temp_segment_final, level_low_gcc, level_high_gcc)

#den
temp_plot = makeplot(data1_1, range, title, metric='density') 
add_window(t_start-50, t_start+50)
color_segments_mat(temp_plot, temp_segment_final, level_low_den, level_high_den)


