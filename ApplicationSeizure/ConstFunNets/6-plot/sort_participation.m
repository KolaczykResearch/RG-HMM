function nodes_ordered = sort_participation(participation, target_com_id)
% sort nodes for nice visualization 

% First sort nodes by length of participation in target community from
% small to large
[~, indices] = sort(sum(participation == target_com_id, 1), 'descend');
participation = participation(:, indices); 

% sort nodes id by the time of participation from later to earlier
[n_row, ~] = size(participation);
nodes_ordered = [];
for row=1:n_row
    % find node ids that participate in at this time
    temp = find(participation(row,:) == target_com_id);
    % append only the non-duplated node id
    nodes_ordered = [nodes_ordered, setdiff(temp, nodes_ordered)];
end
% 
if length(nodes_ordered)<length(indices)
    nodes_ordered = [nodes_ordered, setdiff(1:length(indices), nodes_ordered)];
end
% participation = participation(:, fliplr(nodes_ordered)); 
nodes_ordered = indices(fliplr(nodes_ordered));
end
