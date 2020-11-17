function file = nets_filename(cfg)

params_str = prepdata_filename(cfg);
idx = strfind(params_str, '_');
params_str = [params_str(idx+1:end) '_' cfg.infer.method];
if cfg.infer.scale == "global"
    params_str = [params_str '_scaled'];
elseif cfg.infer.scale == "win"
    params_str = [params_str '_scaledwin'];
end
if cfg.infer.smooth
    params_str = [params_str '_smoothed'];
end
file = ['nets_' params_str];

end