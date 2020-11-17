function cfg = preprocess_data(cfg)
% Preprocess the data

global dynanets_default;

% for i = 1 : length(cfg.data.patients)
%     pat = cfg.data.patients{i};
%     for j = 1 : length(cfg.data.seizures{i})
%         sz = cfg.data.seizures{i}{j};
% 
%     end
% end
%         
% fprintf(['2-preprocess: ' pat '_' sz '\n']);

fprintf(['2-preprocess: ' dynanets_default.inputdatapath '\n']);

% prepdata_file = [dynanets_default.prepdatapath '/' pat '_' sz '/' prepdata_filename(cfg) '.mat'];

prepdata_file = dynanets_default.prepdatapath;

fprintf('... reading in the build file. \n')
% s = load('../../../../dppm/Output_Data/Simulation_Merge_1/data_pad[0,0].mat')
% load in EcoG, Name, Patient, MatlabFile ...
s = load([dynanets_default.inputdatapath]); %load in data = [time, electrodes], hdr

% remove NaN columns from bad channels
s.data = s.data(:, ~any(isnan(s.data), 1));

% subset rows of data to include only 120s before seizure onset and ending 30s after seizure termination
fs = s.hdr.info.sfreq; 
t = cumsum(repmat(1/fs, 1, size(s.data, 1)));
left = max(cfg.infer.szstart -120, 0);
right = cfg.infer.szend + 30;
s.data = s.data(t>left & t<right, :);
s.time = t(t>left & t<right);

% subset columns of data to include only left/right semisphere
if strcmp(cfg.infer.semisphere, 'left') == 1
    fprintf('..semisphere: left\n');
    s.data = s.data(:,1:106);
end
if strcmp(cfg.infer.semisphere, 'right') == 1
    fprintf('..semisphere: right\n');
    s.data = s.data(:,107:176);
end

% Reference data
if strcmp(cfg.preprocess.ref, 'cavg') == 1
    fprintf('... common average reference the data. \n')
    re_referenced_data = commonAvgRef(s.data);
    %diagnostics_car(s.ECoG.Time, s.ECoG.Data, re_referenced_data, [dynanets_default.outfigpath '/' pat '_' sz '/fig/referencing_diagnositics']);
    s.data = re_referenced_data;
    clear re_referenced_data; 
end

% Filter data
% if strcmp(cfg.preprocess.filt, 'firls') == 1
%     fprintf('... filtering the data. \n')
%     fs = s.ECoG.SamplingRate; 
%     filtered_data = lsfilter(s.ECoG.Data, fs, cfg.preprocess.band); % s.ECoG.Data [50000?64 double]
%     %diagnostics_filter(s.ECoG.Time, s.ECoG.Data, filtered_data, fs, cfg.preprocess.band, [dynanets_default.outfigpath '/' pat '_' sz '/fig/filter_diagnositics']);
%     s.ECoG.Data   = filtered_data;
%     clear filtered_data;
% end

if strcmp(cfg.preprocess.filt, 'firls') == 1
    fprintf('... filtering the data. \n')
    fs = s.hdr.info.sfreq; 
    filtered_data = lsfilter(s.data, fs, cfg.preprocess.band); % s.data [times, electrodes]
    %diagnostics_filter(s.ECoG.Time, s.ECoG.Data, filtered_data, fs, cfg.preprocess.band, [dynanets_default.outfigpath '/' pat '_' sz '/fig/filter_diagnositics']);
    s.data  = filtered_data;
    clear filtered_data;
end

%             % Save a simple representation of the whole data for diagnostics
%             f = figure('Visible', 'off');
%             plotchannels(s.ECoG.Time, s.ECoG.Data);
%             print(f, '-dpng', '-r200',  [data_file '.png']);
%             close(f)

% Save preprocessed data in a different file
fprintf('... saving the preprocessed data. \n')
save(prepdata_file, '-struct', 's');
clear s;

end
