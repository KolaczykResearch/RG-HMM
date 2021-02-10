function plot_gcc(nets, szstart, szend)
% Plot the density.

T = size(nets.C,3);             % Total # time points.
n = size(nets.C,1);             % # of nodes
gcc = zeros(1,T);               % Output gcc variable.
for t=1:T                       % For each time index,
    C0 = nets.C(:,:,t);         % ... get the network,
    gcc0 = length(largestcomponent(C0))/n;   % ... compute the gcc, 
    gcc(t) = gcc0;              % ... save it.
end       

plot(nets.t, gcc);
ylim([0,1])
hold on
plot([szstart szstart], ylim, 'r-');
plot([szend szend], ylim, 'r-');
hold off 

axis tight
xlabel('Time (s)')
ylabel('Proportion of nodes in GCC')
end            
