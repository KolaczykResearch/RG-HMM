function file = prepdata_filename(cfg)

params_str = data_filename(cfg);
if strcmp(cfg.preprocess.ref, 'cavg') || strcmp(cfg.preprocess.ref, 'bipolar') == 1
    params_str  = [params_str '_' cfg.preprocess.ref];
end
if strcmp(cfg.preprocess.filt, 'firls') == 1
    params_str = [params_str '_[' num2str(cfg.preprocess.band(1)) '-' num2str(cfg.preprocess.band(2)) ']Hz'];
end

file = ['prep' params_str];


end