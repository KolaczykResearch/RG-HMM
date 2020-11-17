function vnets = network_temporal_vote(nets, win)

% we want odd win size for simplicity
assert(mod(win, 2) ~= 0)

half_width = floor(win / 2);
vnets = nets;

for i = half_width + 1 : length(nets) - half_width
    vnets(:,:,i) = median(nets(:,:,i - half_width : i + half_width), 3);
%     vnets(:,:,i) = median(vnets(:,:,i - half_width : i + half_width), 3);
end

end