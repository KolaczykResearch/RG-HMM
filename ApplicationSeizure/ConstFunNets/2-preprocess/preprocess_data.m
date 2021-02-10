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

% remove NaN columns from bad channels, remove bad channel names
s.hdr.info.ch_names = s.hdr.info.ch_names(~any(isnan(s.data), 1));
s.data = s.data(:, ~any(isnan(s.data), 1));

% if true, subset rows of data to include only 120s before seizure onset and ending 30s after seizure termination
if cfg.infer.subset
    fs = s.hdr.info.sfreq; 
    t = cumsum(repmat(1/fs, 1, size(s.data, 1)));
    left = max(cfg.infer.szstart -120, 0);
    right = cfg.infer.szend + 30;
    s.data = s.data(t>left & t<right, :);
    s.time = t(t>left & t<right);
else
    fs = s.hdr.info.sfreq; 
    t = cumsum(repmat(1/fs, 1, size(s.data, 1)));
    s.time = t;
end

% subset columns of data to include only left/right semisphere
if strcmp(cfg.infer.hemisphere, 'left') == 1
    fprintf('..hemisphere: left\n');
    s.data = s.data(:,1:106);
    s.hdr.info.ch_names = s.hdr.info.ch_names(1:106);
end
if strcmp(cfg.infer.hemisphere, 'right') == 1
    fprintf('..hemisphere: right\n');
    s.data = s.data(:,107:176);
    s.hdr.info.ch_names = s.hdr.info.ch_names(107:176);
end

% Reference data
if strcmp(cfg.preprocess.ref, 'cavg') == 1
    fprintf('... common average reference the data. \n');
    re_referenced_data = commonAvgRef(s.data);
    %diagnostics_car(s.ECoG.Time, s.ECoG.Data, re_referenced_data, [dynanets_default.outfigpath '/' pat '_' sz '/fig/referencing_diagnositics']);
    s.data = re_referenced_data;
    clear re_referenced_data; 
end

% manually code up bipolar reference for left hemisphere or whole brain only!
% this does not apply to right hemesphere
if strcmp(cfg.preprocess.ref, 'bipolar') == 1
    fprintf('... bipolar reference the data. \n');
    re_referenced_data = [];
    %left
    re_referenced_data = [re_referenced_data, s.data(:, 1:15)-s.data(:, 2:16)];
    re_referenced_data = [re_referenced_data, s.data(:, 17:29)-s.data(:, 18:30)];
    re_referenced_data = [re_referenced_data, s.data(:, 31:45)-s.data(:, 32:46)];
    re_referenced_data = [re_referenced_data, s.data(:, 47:57)-s.data(:, 48:58)];
    re_referenced_data = [re_referenced_data, s.data(:, 59:71)-s.data(:, 60:72)];
    re_referenced_data = [re_referenced_data, s.data(:, 73:83)-s.data(:, 74:84)];
    re_referenced_data = [re_referenced_data, s.data(:, 85:93)-s.data(:, 86:94)];
    re_referenced_data = [re_referenced_data, s.data(:, 95:105)-s.data(:, 96:106)];
    %right
    if strcmp(cfg.infer.hemisphere, 'whole') == 1
        re_referenced_data = [re_referenced_data, s.data(:, 107:121)-s.data(:, 108:122)];
        re_referenced_data = [re_referenced_data, s.data(:, 123:135)-s.data(:, 124:136)];
        re_referenced_data = [re_referenced_data, s.data(:, 137:149)-s.data(:, 138:150)];
        re_referenced_data = [re_referenced_data, s.data(:, 151:163)-s.data(:, 152:164)];
        re_referenced_data = [re_referenced_data, s.data(:, 165:175)-s.data(:, 166:176)];
    end
    s.data = re_referenced_data;
    clear re_referenced_data; 
    % get ch_names_new
    % left
    s.hdr.info.ch_names_new = cellstr(strcat(s.hdr.info.ch_names(1:15), '-', string(regexp(string(s.hdr.info.ch_names(2:16)), '\d*', 'match'))));
    temp = cellstr(strcat(s.hdr.info.ch_names(17:29), '-', string(regexp(string(s.hdr.info.ch_names(18:30)), '\d*', 'match'))));
    s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;
    temp = cellstr(strcat(s.hdr.info.ch_names(31:45), '-', string(regexp(string(s.hdr.info.ch_names(32:46)), '\d*', 'match'))));
    s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;
    temp = cellstr(strcat(s.hdr.info.ch_names(47:57), '-', string(regexp(string(s.hdr.info.ch_names(48:58)), '\d*', 'match'))));
    s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;
    temp = cellstr(strcat(s.hdr.info.ch_names(59:71), '-', string(regexp(string(s.hdr.info.ch_names(60:72)), '\d*', 'match'))));
    s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;
    temp = cellstr(strcat(s.hdr.info.ch_names(73:83), '-', string(regexp(string(s.hdr.info.ch_names(74:84)), '\d*', 'match'))));
    s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;
    temp = cellstr(strcat(s.hdr.info.ch_names(85:93), '-', string(regexp(string(s.hdr.info.ch_names(86:94)), '\d*', 'match'))));
    s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;
    temp = cellstr(strcat(s.hdr.info.ch_names(95:105), '-', string(regexp(string(s.hdr.info.ch_names(96:106)), '\d*', 'match'))));
    s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;
    % right
    if strcmp(cfg.infer.hemisphere, 'whole') == 1
        temp = cellstr(strcat(s.hdr.info.ch_names(107:121), '-', string(regexp(string(s.hdr.info.ch_names(108:122)), '\d*', 'match'))));
        s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;
        temp = cellstr(strcat(s.hdr.info.ch_names(123:135), '-', string(regexp(string(s.hdr.info.ch_names(124:136)), '\d*', 'match'))));
        s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;
        temp = cellstr(strcat(s.hdr.info.ch_names(137:149), '-', string(regexp(string(s.hdr.info.ch_names(138:150)), '\d*', 'match'))));
        s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;    
        temp = cellstr(strcat(s.hdr.info.ch_names(151:163), '-', string(regexp(string(s.hdr.info.ch_names(152:164)), '\d*', 'match'))));
        s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;
        temp = cellstr(strcat(s.hdr.info.ch_names(165:175), '-', string(regexp(string(s.hdr.info.ch_names(166:176)), '\d*', 'match'))));
        s.hdr.info.ch_names_new(length(s.hdr.info.ch_names_new)+(1:length(temp))) = temp;
    end
    %exclude channels 
    ch_exclude = {'LOF1-2'; 'LOF2-3'; 'LOF3-4'; 'LOF4-5'; 'LOF5-6'; 'LOF6-7'; 'LOF7-8'; 'LOF8-9'; 'LOF9-10'; 'LMF11-12'; 'LMF12-13'; 'LMF13-14'; 'LI11-12'; 'LMT11-12'; 'LMT12-13'; 'LMT13-14'; 'LPT7-8'; 'LPT8-9';'LPT9-10';'LPT10-11';'LPT11-12'; 'LAL4-5'; 'LAL9-10'; 'LAL10-11';'LAL11-12'};
    exclude_bool = 0;
    for i=1:length(ch_exclude)
        ch = ch_exclude(i);
        exclude_bool = exclude_bool + strcmp(s.hdr.info.ch_names_new, ch);
    end
    s.data = s.data(:,~logical(exclude_bool));
    s.hdr.info.ch_names_new = s.hdr.info.ch_names_new(~logical(exclude_bool));
end

% Filter data
if strcmp(cfg.preprocess.filt, 'firls') == 1
    fprintf('... filtering the data. \n')
    fs = s.hdr.info.sfreq; 
    filtered_data = lsfilter(s.data, fs, cfg.preprocess.band); % s.data [times, electrodes]
    %diagnostics_filter(s.ECoG.Time, s.ECoG.Data, filtered_data, fs, cfg.preprocess.band, [dynanets_default.outfigpath '/' pat '_' sz '/fig/filter_diagnositics']);
    s.data  = filtered_data;
    clear filtered_data;
end

% Save preprocessed data in a different file
fprintf('... saving the preprocessed data. \n')
save(prepdata_file, '-struct', 's');
clear s;

end
