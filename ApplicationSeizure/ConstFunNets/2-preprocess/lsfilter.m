function  y = lsfilter(data,fs,BAND)
% LSFILTER creates FIRLS filter using least-squares error minimization. 
% Filtfilt is applied to data so that there is no phase shift.
% Data should be NxM with N the number time step and M the number of channels. BAND is the
% band of frequencies to be included in the data, all other frequencies
% will be excluded.  Fs is the sampling frequency. Y is filtered data.

ORDER = 300; % needs 3*order samples
order = min(ORDER, floor(size(data,1)/3));
if order ~= ORDER
    warning(['lsfilter: using filter order ' num2str(order) ' instead of ' num2str(ORDER)]);
end

fNQ = 0.5 * fs;             % Nyquist frequency, half the sampling rate

if BAND(1) > 0
    if BAND(2) < fNQ
        f = [0 (BAND(1)-1)/fNQ BAND/fNQ (BAND(2)+1)/fNQ 1];
        a = [0 0 1 1 0 0];
    elseif BAND(2) == fNQ
        f = [0 (BAND(1)-1)/fNQ BAND/fNQ];
        a = [0 0 1 1];
    else
        f = [0 (BAND(1)-1)/fNQ BAND(1)/fNQ 1];
        a = [0 0 1 1];
    end
elseif BAND(1) == 0
    if BAND(2) < fNQ
        f = [BAND/fNQ (BAND(2)+1)/fNQ 1];
        a = [1 1 0 0];
    elseif BAND(2) == fNQ
        f = [BAND/fNQ];
        a = [1 1];
    else
        f = [BAND(1)/fNQ 1];
        a = [1 1];
    end
        
end

flt = firls(order,f,a);
y = filtfilt(flt,1,data);

end

