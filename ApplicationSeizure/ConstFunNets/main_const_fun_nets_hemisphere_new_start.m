cfg = [];

% % Build simulation data
% cfg.data.run = true;
% cfg.data.patients = {'Simulation'};
% cfg.data.simscenario = 'Merge';         % Choice: 'Null', 'Expand', 'Contract', 'Split', 'Merge'
% cfg.data.nbsim = 1;                     % Run a single simulation of one "seizure".
% cfg.data.seizures = {{[cfg.data.simscenario '_1']}};
% cfg.data.padding = [0 0];

% Preprocessing settings
cfg.preprocess.run = true;
cfg.preprocess.ref = '';                % Choice: '', 'cavg', 'bipolar'
cfg.preprocess.filt = 'firls';          % Linear-phase FIR filter
cfg.preprocess.band = [4 50];           % Bandpass

% Inference settings
cfg.infer.run = true;
cfg.infer.method = 'corr_0_lag';        % Choice: 'corr', 'corr_0_lag'
cfg.infer.windowsize = 1;               % Window size for net inference (s). 1, 0.5
cfg.infer.windowstep = 0.5;             % Window overlap (s), 0.5, 0.2, 0.1
cfg.infer.smooth = false;               % Use vote to smooth networks in time
%cfg.infer.scale  = true;                % Scale variance of correlation using all time.
cfg.infer.scale  = 'global';              % Scale variance of correlation using all time.

cfg.infer.semisphere = 'left';
%cfg.infer.semisphere = 'right';


% cfg.infer.szstart = 64;                 % seizure begins near X second. 
% cfg.infer.szend = 314;                   % seizure ends near Y second.
% cfg.infer.id = 'P1S1_hemi';   

% cfg.infer.szstart = 203;                 % seizure begins near X second. 
% cfg.infer.szend = 412;                   % seizure ends near Y second.
% cfg.infer.id = 'P1S2_hemi';   

cfg.infer.szstart = 80;                  % seizure begins near X second. 
cfg.infer.szend = 250;                   % seizure ends near Y second.
cfg.infer.id = 'P1S3_hemi';  

% infer networks beginning 120s before seizure onset and ending 30s after seizure termination

global dynanets_default;
% dynanets_default.inputdatapath = 'DataSeizure/EP001_clip_1_ecog.mat';
% dynanets_default.inputdatapath = 'DataSeizure/EP001_clip_2_ecog.mat';
dynanets_default.inputdatapath = 'DataSeizure/EP001_clip_3_ecog.mat';

[filepath, name, ext] = fileparts(dynanets_default.inputdatapath);

dynanets_default.prepdatapath = ['DataSeizure/prepdata/' name '_ws_' num2str(cfg.infer.windowsize) '_wstp_' num2str(cfg.infer.windowstep) '_sub_hemi_' cfg.infer.semisphere '.mat'];
dynanets_default.netsdatapath = ['Data/nets_' name '_ws_' num2str(cfg.infer.windowsize) '_wstp_' num2str(cfg.infer.windowstep) '_sub_hemi_' cfg.infer.semisphere '.mat']; % output

dynanets_default.codepath = mfilename('fullpath');
dynanets_default.codepath = fileparts(dynanets_default.codepath);
addpath([dynanets_default.codepath '/2-preprocess']);
addpath([dynanets_default.codepath '/3-infer']);

cfg = main_dynanets(cfg);
