# ========================================================
# == Automatic Segment Finder ============================
# ========================================================
# Given a ROI (identified from the large dynamic community), return the 
# left and right indices for all ramp-ups in this region using one metric.
# == Input:
# Y: network time series, list of matrices
# t_obs: sequence of times
# t_start: seizure onset
# ROI: 2 dimensional vector 
# metric: 'gcc' or 'density'
# == Output:
# segment(s) of ramp-up(s)
# segments_all: vector, left and right indices for each low-high-low segment
# segments_mat: matrix, each row for each low-high segment
auto_segment_finder = function(Y, t_obs, t_start, ROI, metric='gcc'){
  N = nrow(Y[[1]])
  if(is.null(ROI)){
    # 50s before t_start, and 50s after t_start
    index_left = get_index(t_start-50, t_obs)
    index_right = get_index(t_start+50, t_obs)
  }else{
    if(is.vector(ROI)){
      index_left = get_index(ROI[1], t_obs)
      index_right = get_index(ROI[2], t_obs)
    }else{
      stop('Only accept one ROI.')
    }
  }
  # get gcc or den
  gcc_or_den = numeric(0)
  for (t in index_left:index_right){
    if(metric=='gcc'){
      gcc_or_den = c(gcc_or_den, getSizeGCC(Y[[t]], order=1)/N)
    }else{
      gcc_or_den = c(gcc_or_den, get_density(Y[[t]], N))
    }
  }
  # split a curve into different curve segments by 
  # finding low_high_low segments of gcc_or_den
  # --> set hard limit to q1 and q3
  if(metric=='gcc'){
    level_low = min(0.14, quantile(gcc_or_den, 0.25))
    # level_low = min(0.15, quantile(gcc_or_den, 0.25))
    level_high = max(0.8, quantile(gcc_or_den, 0.75))
  }else{
    # level_low = min(0.05, quantile(gcc_or_den, 0.25))
    level_low = min(0.06, quantile(gcc_or_den, 0.25))
    level_high = max(0.15, quantile(gcc_or_den, 0.75))
  }
  level_high = min(0.7*max(gcc_or_den), level_high)
  level_low = max(1.1*min(gcc_or_den), level_low)
  # level_low = max(2*min(gcc_or_den), level_low)
  # --> find low-high-low
  set = find_low_high_low_segments_splits(gcc_or_den, level_low, level_high)
  set = set + index_left - 1 
  
  # for each segment, trim off the high-low tail 
  segments_mat = NULL
  for(i in 2:length(set)){
    a = set[i-1]
    b = set[i]
    while(gcc_or_den[b-index_left+1]<level_high){
      b = b-1
    }
    if((b-a+1)>=8){
      segments_mat = rbind(segments_mat, c(a,b)) # only include segments with duration > 4s
    }
  }
  return(list(segments_all=set, segments_mat=segments_mat, 
              level_low=level_low, level_high=level_high))
}


# ========================================================
# == Automatic Segment Finder with Fine Tuning ===========
# ========================================================
# Given a ROI (identified from the large dynamic community), return the 
# left and right indices for all fine-tuned ramp-ups in this region using two metrics.
# == Input:
# Y: network time series, list of matrices
# t_obs: sequence of times
# t_start: seizure onset
# ROI: 2 dimensional vector 
# metric: 'gcc' or 'density'
# == Output:
# segment(s) of ramp-up(s)
# segments_all: matrix, each row for each low-high segment
# segments_fine_tuned: matrix, each row for each fine-tuned low-high segment
# level_low_gcc, level_high_gcc, level_low_den, level_high_den
auto_segment_finder_fine_tune = function(Y, t_obs, t_start=NULL, ROI=NULL){
  if(is.null(t_start) & is.null(ROI)){stop('Either t_start or ROI is needed!')}
  N = nrow(Y[[1]])
  gcc = sapply(Y, getSizeGCC, order=1)/N
  den = sapply(Y, get_density, N=N)
  # find low-high segments with gcc 
  temp_segment = auto_segment_finder(Y, t_obs, t_start, ROI, metric='gcc')
  segments_mat_gcc = temp_segment$segments_mat
  level_low_gcc = temp_segment$level_low
  level_high_gcc = temp_segment$level_high
  # find low-high segments with density 
  temp_segment = auto_segment_finder(Y, t_obs, t_start, ROI, metric='density')
  segments_mat_den = temp_segment$segments_mat
  level_low_den = temp_segment$level_low
  level_high_den = temp_segment$level_high
  # find overlapping
  segments_mat_overlapped = find_overlapping(segments_mat_gcc, segments_mat_den)
  # fine tune
  # --> remove segments that are < 4s
  bools = (segments_mat_overlapped[,2]-segments_mat_overlapped[,1]+1)>=8
  segments_mat_overlapped = segments_mat_overlapped[bools,]
  if(!is.matrix(segments_mat_overlapped)) segments_mat_overlapped = matrix(segments_mat_overlapped, nrow=1, ncol=2)
  # --> exclude those with high start and low end
  segments_mat_final = NULL
  for(i in 1:nrow(segments_mat_overlapped)){
    left = segments_mat_overlapped[i,1]
    right = segments_mat_overlapped[i,2]
    cond1 = gcc[left]<=level_low_gcc
    cond2 = gcc[right]>=0.5
    cond3 = den[left]<=(level_low_den*1.6)
    if (cond1 & cond2 & cond3){
      segments_mat_final = rbind(segments_mat_final, segments_mat_overlapped[i,])
    }
  }
  return(list(segments_fine_tuned=segments_mat_final, segments_all=segments_mat_overlapped,
              level_low_gcc=level_low_gcc, level_high_gcc=level_high_gcc,
              level_low_den=level_low_den, level_high_den=level_high_den))
}

# ========================================================
# == Define helper functions =============================
# ========================================================
# get the index in the sequence of time t in seconds
get_index = function(t, t_all){
  which(t_all >= t)[1]
}

# get the size of gcc and 1st/2nd connected component 
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
      return (sort(clusters_size)[id])
    }else{
      return (0)
    }
  }else{
    cat ('Order Not Defined!')
  }
  return(0)
}

# get density
get_density = function(net_t, N){
  nedge = sum(net_t==1)/2
  den = nedge/(N*(N-1)/2)
  return(den)
}

# function called in auto_segment_finder()
# Divide a given time series into low-high-low segments
find_low_high_low_segments_splits = function(gcc_or_den, level_low, level_high){
  set = NULL
  # Anchor the left point of a potential segment as the first point below level_low
  anchor = 1
  while(gcc_or_den[anchor]>level_low){
    anchor = anchor + 1
  }
  set = c(set, anchor)
  # Start segmenting time series with sliding window
  while(set[length(set)] < length(gcc_or_den)){
    i = 1
    while((gcc_or_den[anchor+i]<level_high) & (anchor+i<length(gcc_or_den))){
      i = i+1 # grow the segment until it reaches Q3
    }
    while((gcc_or_den[anchor+i]>level_low) & (anchor+i<length(gcc_or_den))){
      i = i+1 # keep growing the segment until it reaches Q1
    }
    anchor = anchor + i
    set = c(set, anchor)
  }
  return(set)
}

# Given two sets of segments, return the overlapping 
find_overlapping = function(segments_mat_gcc, segments_mat_den){
  # matrix: 1st col: splitting points; 2nd col: left(1) or right(2) point of segment; 3rd col: from gcc(1) or den(2)
  all_points = matrix(0, (nrow(segments_mat_gcc)+nrow(segments_mat_den))*2, 3) 
  left_end =max(segments_mat_gcc[1,1], segments_mat_den[1,1])
  right_end = min(segments_mat_gcc[nrow(segments_mat_gcc),2], segments_mat_den[nrow(segments_mat_den),2])
  id = 1
  for(i in 1:nrow(segments_mat_gcc)){
    all_points[id,]=c(segments_mat_gcc[i,1], 1, 1)
    all_points[id+1,]=c(segments_mat_gcc[i,2], 2, 1)
    id = id+2
  }
  for(i in 1:nrow(segments_mat_den)){
    all_points[id,]=c(segments_mat_den[i,1], 1, 2)
    all_points[id+1,]=c(segments_mat_den[i,2], 2, 2)
    id = id+2
  }
  all_points = all_points[((all_points[,1]<=right_end) & (all_points[,1]>=left_end)),]
  # sort all splitting points
  all_points = all_points[order(all_points[,1]),]
  # for each consecutive pair, decide if it is overlapped segment
  segments_mat_overlapped = NULL
  for(i in 2:nrow(all_points)){
    if((all_points[i-1, 2]==1) & (all_points[i, 2]==2)){ # one left split, one right
      segments_mat_overlapped = rbind(segments_mat_overlapped, c(all_points[i-1, 1], all_points[i, 1]))
    }
  }
  return(segments_mat_overlapped)
}

