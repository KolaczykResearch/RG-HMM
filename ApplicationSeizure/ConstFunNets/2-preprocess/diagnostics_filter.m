function diagnostics_filter(time, x, y, fs, band, outpath)

% Do some diagnostics comparing x and the filter result y, put in figs folder                                                    

h(1) = subplot(3,1,1);
plot(time,x(:,4),time,y(:,4),'r');                      % Plots signals filtered and unfiltered
xlabel('Time (s)')
legend('Original Signal','Filtered Signal');
title( ['Channel 4 Signal Filtered and Unfiltered [' num2str(band) ' Hz]'])
axis tight

% p = mean(x,1);                                         
% w = p'* ones(1,size(x,1));
% x = x - w';

x = bsxfun(@minus, x, mean(x,1));  % subtract mean from each electrode before spectrum estimation

window = [2 0.5];                                       % Compute the spectrogram.
TW = 2*2;                                               % Time bandwidth product
nbtaper = 2*TW-1;                                       % no. tapers estimates of frequency in window
params = struct('tapers', [TW, nbtaper],'pad',-1,'Fs',fs,'fpass',[0 100]);

[S, t, F] = mtspecgramc(x(:,4),window,params);
t_real = t + time(1) + window(1)/2;
P = 10 * log10(abs(S'));                                            
h(2) = subplot(3,1,2);
surf(t_real,F,P,'EdgeColor','none');                         % Plots log spectrum of unfiltered data
axis xy; colormap(jet); view(0,90); axis tight
mx = max(P(:));  mn = min(P(:));
caxis([mn mx])
title('Unfiltered data spectrogram, channel = 4')

[S, ~, F] = mtspecgramc(y(:,4),window,params);
P = 10 * log10(abs(S'));
h(3) = subplot(3,1,3);
surf(t_real,F,P,'EdgeColor','none');                         % Plots log spectrum of unfiltered data
axis xy; colormap(jet); view(0,90); axis tight
caxis([mn mx])
title('Filtered data spectrogram, channel = 4')

h = [h(1) h(2) h(3)];
linkaxes(h,'x')
print(gcf, '-djpeg', outpath);
% saveas(h, outpath,'fig');
close

end