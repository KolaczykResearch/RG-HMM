cfg = [];

% Preprocessing settings
cfg.preprocess.run = true;
cfg.preprocess.ref = 'bipolar';         % Choice: '', 'cavg', 'bipolar'
cfg.preprocess.filt = 'firls';          % Linear-phase FIR filter
cfg.preprocess.band = [4 50];           % Bandpass

% Inference settings
cfg.infer.run = true;
cfg.infer.method = 'corr_0_lag';        % Choice: 'corr', 'corr_0_lag'
cfg.infer.windowsize = 1;               % Window size for net inference (s). 1, 0.5
cfg.infer.windowstep = 0.5;             % Window overlap (s), 0.5, 0.2, 0.1
cfg.infer.smooth = false;               % Use vote to smooth networks in time
cfg.infer.scale = 'global';             % Scale variance of correlation using all time, choice 'global', 'local'.
cfg.infer.thr = 'dynamic';              % Choice: 'dynamic', 'static'. Convert weighted networks into binary networks with varying/static threshold
cfg.infer.fdr = 0.5;                    % False discovery rate 

cfg.infer.hemisphere = 'left';
% cfg.infer.hemisphere = 'whole';
cfg.infer.subset = true;                % if true, subset rows of data to include only 120s before seizure onset and ending 30s after seizure termination

cfg.infer.szstart = 64;               % seizure begins near X second. 
cfg.infer.szend = 314;                % seizure ends near Y second.
cfg.infer.id = 'P1S1_hemi';   

% cfg.infer.szstart = 203;              % seizure begins near X second. 
% cfg.infer.szend = 412;                % seizure ends near Y second.
% cfg.infer.id = 'P1S2_hemi';   

% cfg.infer.szstart = 80;               % seizure begins near X second. 
% cfg.infer.szend = 250;                % seizure ends near Y second.
% cfg.infer.id = 'P1S3_hemi';  

% REPLACE WITH YOUR PATH FOR inputdatapath, prepdatapath, netsdatapath
global dynanets_default;
dynanets_default.inputdatapath = '../../DataSeizure/EP001_clip_1_ecog.mat';
% dynanets_default.inputdatapath = '../../DataSeizure/EP001_clip_2_ecog.mat';
% dynanets_default.inputdatapath = '../../DataSeizure/EP001_clip_3_ecog.mat';

[filepath, name, ext] = fileparts(dynanets_default.inputdatapath);

dynanets_default.prepdatapath = ['../../DataSeizure/prepdata_fdr/' name '_ws_' num2str(cfg.infer.windowsize) '_wstp_' num2str(cfg.infer.windowstep) '_sub_hemi_' cfg.infer.hemisphere '_' cfg.preprocess.ref '_fdr_' num2str(cfg.infer.fdr) '.mat'];
dynanets_default.netsdatapath = ['../../DataSeizure/networks_fdr/nets_' name '_ws_' num2str(cfg.infer.windowsize) '_wstp_' num2str(cfg.infer.windowstep) '_sub_hemi_' cfg.infer.hemisphere '_' cfg.preprocess.ref '_fdr_' num2str(cfg.infer.fdr) '.mat'];

% REPLACE WITH YOUR PATH
dynanets_default.codepath = 'RG-HMM/ApplicationSeizure/ConstFunNets_bipolar_fdr';
addpath([dynanets_default.codepath '/2-preprocess']);
addpath([dynanets_default.codepath '/3-infer']);

cfg = main_dynanets(cfg);









