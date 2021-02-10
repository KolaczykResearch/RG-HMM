function plot_voltage_trace_hemi(s, szstart, szend)
% Plot the preprocessed voltage traces (node subset) for left hemisphere.

s.data = s.data./max(s.data, [], 1);
% s.data = s.data./max(s.data, [], [1, 2]);

% plot node subset
to_plot = (1:10:size(s.data,2));
for k=to_plot
    plot(s.time, 10*s.data(:,k)+k, 'LineWidth', 0.5, 'Color', 'black') 
    hold on
end
plot([szstart szstart], ylim, 'r-');
plot([szend szend], ylim, 'r-');
hold off
set(gca, 'YTick', to_plot)
set(gca, 'YTickLabel', s.hdr.info.ch_names_new(to_plot))

axis tight
xlabel('Time (s)')
ylabel('')
end

