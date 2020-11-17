% Compute the standard deviation of the cross correlation for the entire
% data set, then average over all time (and electrode pairs) to compute a
% "scale" to use in the network inference.

function [scale] = compute_correlation_scale(d,t,cfg)

  % Divide the data into windows, with overlap.
  i_total = 1+floor((t(end)-t(1)-cfg.infer.windowsize) / cfg.infer.windowstep);       % # intervals.

  % Indices for non-redundant edges in network.
  iupper = find(triu(ones(size(d,2)),1))'; 
  
  % Output variables.
  stdsij   = zeros(i_total,length(iupper));
  
  parfor k=1:i_total                                             %For each window,

      t_start = t(1) + (k-1) * cfg.infer.windowstep; %#ok<PFBNS> %... get window start time [s],
      t_stop  = t_start + cfg.infer.windowsize;                  %... get window stop time [s],
      indices = t >= t_start & t < t_stop;                       %... find indices for window in t,
      
      stdsij0     = compute_correlation_establish_variance0(d(indices,:));
      stdsij(k,:) = stdsij0(iupper);
      
  end
  
  if cfg.infer.scale == "global"
    scale = nanmean(stdsij(:));
  elseif cfg.infer.scale == "win"
    scale = nanmean(stdsij,2);
  end

end

%%%%%%%%
%NOTES.
%
%Compute standard deviation of cross correlation and return.

function [stdsij] = compute_correlation_establish_variance0(dgrid)

  %Define some useful variables.
  winSize = size(dgrid,1);          %The length of the data in time.
  N = size(dgrid,2);                %The number of electrodes.
  
  %Consider max correlations only within [-winSize/4, winSize/4].
  maxLags=floor(winSize/4);
  
  %Do the correlation analysis.    
  X = dgrid(:,:)';                                      %Call the data X.
  X = X - mean(X,2) * ones(1,winSize);                  %Remove the mean.
  X = X ./ (std(X,[],2) * ones(1,winSize));             %Set stdev to 1.
  for k=1:N                                             %Detrend activity at each electrode.
      X(k,:) = detrend(X(k,:));
  end
  
  %Zero pad.
  pad = zeros(N, maxLags);      %The padding has length of maxLags.
  X = [pad, X, pad];            %Add padding to beginning and end of data.
  winSize = size(X,2);          %The new length of the data in time.
  
  ac = ifft( fft(X,[],2) .* conj(fft(X,[],2)), [],2);   %Compute the auto-correlation
  ac = ac(:,1) * ac(:,1)';                              %... at zero lag
  ac = repmat(ac, [1,1,winSize]);                       %... make a 3-dim matrix
  ac = permute(ac, [1,3,2]);                            %... and match dims of cross-correlation.
  
                                                        %Compute the cross-correlation
  X3 = repmat(fft(X,[],2),[1,1,N]);                     %... make 3-dim [N,winSize,N]
  cc = ifft(permute(X3,[3,2,1]).*conj(X3), [], 2);      %... use FTs
  cc = cc ./ sqrt(ac);                                  %... scale by AC.

  sij = 0.5*log((1+cc)./(1-cc));                        %Fisher transform the cc.
  stdsij = squeeze(std(  sij(:,[(1:maxLags+1), ...      %Compute std over Fisher cc's.
                                (end-maxLags+1:end)],:), [], 2));
                            
end
  