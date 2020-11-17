function b = bipolarRefSimulation(data)

  n = size(data,2);
  T = size(data,1);
  b = zeros(T,n-1);
  for k=1:n-1
      %dist0 = sqrt(sum((xy-repmat(xy(k,:),[n,1])).^2,2));
      %dist0(k) = NaN;
      %[~,i0] = min(dist0);
      b(:,k) = data(:,k) - data(:,k+1);
  end

end

