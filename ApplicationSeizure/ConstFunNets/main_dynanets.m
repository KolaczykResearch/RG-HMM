function cfg = main_dynanets(cfg)

% Preprocess data and save to .mat
cfg = preprocess_data(cfg);

% Infer functional networks
cfg = infer_nets(cfg);

end
