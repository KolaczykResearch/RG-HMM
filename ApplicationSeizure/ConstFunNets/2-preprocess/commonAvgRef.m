function b = commonAvgRef( data )
% COMMONAVGREF calculates the common average reference of DATA. DATA should
% be NxM with N the number time step and M the number of channels. B is
% referenced data.
% 

channel_avg = mean(data,2); % Average the data over channels, for each moment in time.
CAR = channel_avg * ones(1,size(data,2)); % Compute the CAR for each electrode, and each moment in time. 
b = data - CAR;             % Subtract the CAR from the data.

% possibly faster
% b = bsxfun(@minus, data, mean(data,2));

end

