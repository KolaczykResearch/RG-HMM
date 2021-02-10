function plot_voltage_trace_hemi_sorted(s, szstart, szend, nodes_ordered)
% Plot the preprocessed voltage traces (node subset) for left hemisphere.
s.data = s.data./max(s.data, [], 1);
nodes_ordered = fliplr(nodes_ordered);
[~, order] = sort(nodes_ordered);
% plot node subset
to_plot = (1:10:size(s.data,2));
[~, order2] = sort(order(to_plot));
to_plot_id = to_plot(order2);
for i=[1:length(to_plot)]
    k = to_plot(i);
    id = to_plot_id(i);
    plot(s.time, 10*s.data(:,id)+k, 'LineWidth', 0.5, 'Color', 'black') 
    hold on
end
plot([szstart szstart], ylim, 'r-');
plot([szend szend], ylim, 'r-');
hold off
set(gca, 'YTick', to_plot)
set(gca, 'YTickLabel', s.hdr.info.ch_names_new(to_plot_id))

axis tight
xlabel('Time (s)')
ylabel('')
end


