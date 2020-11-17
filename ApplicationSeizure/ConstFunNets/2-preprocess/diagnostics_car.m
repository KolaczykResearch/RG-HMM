function diagnostics_car(time, x, y, outpath)

% Do some diagnostics comparing x and the re-reference result y, put in figs folder

g = figure();
plot(time,x(:,4),time,y(:,4),'r');                                  % Plots signals before and after re-referencing
xlabel('Time (s)')
legend('Original Signal','Re-referenced Signal');
title('Channel 4 Signal')
axis tight

print(g, '-djpeg', outpath);
% saveas(g, outpath, 'fig');
close;

end