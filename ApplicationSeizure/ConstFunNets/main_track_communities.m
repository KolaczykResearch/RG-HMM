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
fig_name = 'seizure1_left_fdr.pdf';
% fig_name = 'seizure2_left_fdr.pdf';
% fig_name = 'seizure3_left_fdr.pdf';

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
% Plot voltage trace 
fprintf('Plot voltage traces...\n');
voltage = load(Seizure_voltage_path);
subplot(4,1,1)
plot_voltage_trace_hemi(voltage, s.cfg.infer.szstart, s.cfg.infer.szend);

% % Plot the number of coms through time
% fprintf('Plot number of coms...\n');
% subplot(4,1,2)
% plot(t, stats.nb_com);
% ylim([0, max(stats.nb_com)+1]);
% xlabel('Time (s)')
% ylabel('Number of Communities')

% Community participation over time.
fprintf('Plot com participation...\n');
subplot(4,1,2)
cc_nowhite = colorcube;
cc_nowhite = cc_nowhite(1:end-1,:);
participation = stats.participation;
colormap(cc_nowhite)
imagescwithpcolor(t, (1:size(C,1)), participation')
xlabel('Time (s)')
ylabel('Node')

% Plot density over time
fprintf('Plot density...\n');
subplot(4,1,3)
plot_density(nets);

% Plot gcc over time
fprintf('Plot gcc...\n');
subplot(4,1,4)
plot_gcc(nets, s.cfg.infer.szstart, s.cfg.infer.szend);

orient(gcf,'landscape')
print(gcf, fig_name, '-dpdf', '-r1000')


