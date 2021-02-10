% Run dppm with default settings given a time-indexed binary network.
%
% Given C   = [nodes, nodes, time] = time-indexed binary network.
%       t   = [1, time]            = sequence of time
%       k,m = 2,4                  = for a 2-plex of size 4. 

clear

%---- Load toolboxes ------------------------------------------------------

% REPLACE WITH YOUR PATH.
localpath = '../../';
BCT_toolbox_path = [localpath 'dppm_root_dir/BCT/'];
DPP_toolbox_path = [localpath 'dppm_root_dir/dynamic-plex-propagation/'];
DPPM_toolbox_path= [localpath 'dppm_root_dir/dppm/'];

% ------------------------------------------
% PATH for one of seizures 1, 2, 3: binary network
% ------------------------------------------
Seizure_nets_path = [localpath 'DataSeizure/networks_fdr/nets_EP001_clip_1_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat'];
% Seizure_nets_path = [localpath 'DataSeizure/networks_fdr/nets_EP001_clip_2_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat'];
% Seizure_nets_path = [localpath 'DataSeizure/networks_fdr/nets_EP001_clip_3_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat'];

% ------------------------------------------
% PATH for one of seizures 1, 2, 3: preprocessed voltage traces
% ------------------------------------------
Seizure_voltage_path = [localpath 'DataSeizure/prepdata_fdr/EP001_clip_1_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat'];
% Seizure_voltage_path = [localpath 'DataSeizure/prepdata_fdr/EP001_clip_2_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat'];
% Seizure_voltage_path = [localpath 'DataSeizure/prepdata_fdr/EP001_clip_3_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat'];

% ------------------------------------------
% Output directory
% ------------------------------------------
fig_name = 'seizure1_left_fdr_color_sorted.pdf';
% fig_name = 'seizure2_left_fdr_color_sorted.pdf';
% fig_name = 'seizure3_left_fdr_color_sorted.pdf';

% ------------------------------------------
% Path to dependencies
% ------------------------------------------
% Addpath to dpp                            (https://github.com/nathanntg/dynamic-plex-propagation)
addpath(genpath(DPP_toolbox_path));
% Addpath to Brain Connectivity Toolbox     (https://sites.google.com/site/bctnet/)
addpath(genpath(BCT_toolbox_path));
% Addpath to compute community statistics   (https://github.com/Eden-Kramer-Lab/dppm)
addpath(genpath([DPPM_toolbox_path '5-analyze/']))
% Addpath to compute gcc 
addpath('6-plot/')

%---- Load & format example data ------------------------------------------

% REPLACE WITH YOUR DATA.
% load([DPPM_toolbox_path 'simplified_example/C_example.mat']);                 % Load data.
s = load(Seizure_nets_path);
t = transpose(s.nets.t);
C = s.nets.C;

nets = [];  nets.C = C;  nets.t = t;        % Format it for DPPM.

%---- Run DPPM ------------------------------------------------------------
k = 2; m = 4;                               % Set default DPPM parameters
[track.vertices, track.communities] = dpp(C, k,m);

%---- Compute community statistics ----------------------------------------
stats = community_stats(track);

%---- Plot the results ----------------------------------------------------
% % ------------------------
% % Plot unsorted voltage trace 
% % ------------------------
% fprintf('Plot voltage traces...\n');
% voltage = load(Seizure_voltage_path);
% subplot(4,1,1)
% plot_voltage_trace_hemi(voltage, s.cfg.infer.szstart, s.cfg.infer.szend);

% ------------------------
% Community participation over time.
% ------------------------
fprintf('Plot com participation...\n');
subplot(4,1,2)
cc_nowhite = colorcube;
cc_nowhite = cc_nowhite(1:end-1,:);
participation = stats.participation;
cc_nowhite = cc_nowhite(randperm(size(cc_nowhite, 1)), :);

% Sort node ids based on the order of recruitment to the target community
% provide target_com_id for each seizure
% [~, ids] = sort(stats.com_cum_size)
% Seizure 1
target_com_id = 68;
% Seizure 2
% target_com_id = 164;
% Seizure 3
% target_com_id = 116;

nodes_ordered = sort_participation(participation, target_com_id);
participation = participation(:, nodes_ordered);
colormap(cc_nowhite)
imagescwithpcolor(t, (1:size(C,1)), participation')
xlabel('Time (s)')
ylabel('Node')

% ------------------------
% Plot sorted voltage trace 
% ------------------------
fprintf('Plot voltage traces...\n');
voltage = load(Seizure_voltage_path);
subplot(4,1,1)
plot_voltage_trace_hemi_sorted(voltage, s.cfg.infer.szstart, s.cfg.infer.szend, nodes_ordered);

% ------------------------
% Plot density over time
% ------------------------
fprintf('Plot density...\n');
subplot(4,1,3)
plot_density(nets);

% ------------------------
% Plot gcc over time
% ------------------------
fprintf('Plot gcc...\n');
subplot(4,1,4)
plot_gcc(nets, s.cfg.infer.szstart, s.cfg.infer.szend);


% Replot gcc plot subplot(4,1,4) for bold curve
subplot(4,1,4)
T = size(nets.C,3);                          % Total # time points.
n = size(nets.C,1);                          % # of nodes
gcc = zeros(1,T);                            % Output gcc variable.
for t=1:T                                    % For each time index,
    C0 = nets.C(:,:,t);                      % ... get the network,
    gcc0 = length(largestcomponent(C0))/n;   % ... compute the gcc, 
    gcc(t) = gcc0;                           % ... save it.
end       
plot(nets.t, gcc);
ylim([0,1])
hold on
plot([s.cfg.infer.szstart s.cfg.infer.szstart], ylim, 'r-');
plot([s.cfg.infer.szend s.cfg.infer.szend], ylim, 'r-');

% Seizure 1
plot(nets.t(126:134), gcc(126:134), 'LineWidth', 1.5, 'Color', [54, 55, 150]/255);
plot(nets.t(453:463), gcc(453:463), 'LineWidth', 1.5, 'Color', [54, 55, 150]/255);
plot(nets.t(609:624), gcc(609:624), 'LineWidth', 1.5, 'Color', [54, 55, 150]/255);

% Seizure 2
% plot(nets.t(236:253), gcc(236:253), 'LineWidth', 1.5, 'Color', [54, 55, 150]/255);
% plot(nets.t(295:318), gcc(295:318), 'LineWidth', 1.5, 'Color', [54, 55, 150]/255);
% plot(nets.t(627:646), gcc(627:646), 'LineWidth', 1.5, 'Color', [54, 55, 150]/255);

% Seizure 3
% plot(nets.t(229:238), gcc(229:238), 'LineWidth', 1.5, 'Color', [54, 55, 150]/255);
% plot(nets.t(256:269), gcc(256:269), 'LineWidth', 1.5, 'Color', [54, 55, 150]/255);
% plot(nets.t(309:316), gcc(309:316), 'LineWidth', 1.5, 'Color', [54, 55, 150]/255);
% plot(nets.t(385:395), gcc(385:395), 'LineWidth', 1.5, 'Color', [54, 55, 150]/255);

hold off
axis tight
xlabel('Time (s)')
ylabel('Proportion of nodes in GCC')

% Save figure
orient(gcf,'landscape')
print(gcf, fig_name, '-dpdf', '-r1000')
