function ieeg_2_sync_data_recon(sid)
% IEEG_2_SYNC_DATA_RECON - synchronize raw iEEG data and electrode 
% reconstruction. Match channel labels and save channel x-y-z coordinates 
% in data structure.
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Inputs:
% sid = subject ID (e.g., 'NM01')
%
% Example:
% ieeg_2_sync_data_recon('NM01')
%
% Copyright (c) 2022-2025
% EL Johnson, PhD

clearvars -except sid

% set directories
pth = fullfile(pwd, sid);
datdir_data = fullfile(pth, 'Preprocessing'); % output of ieeg_0_read_data
tmp = dir(fullfile(pth, 'Recon', 'Recon_*'));
datdir_recon = fullfile(pth, 'recon', tmp.name, 'FT_Pipeline', ...
    'Electrodes'); % ouput of ieeg_1_recon
clear tmp

savdir = datdir_data;

% load data
data = load(fullfile(datdir_data, [sid '_data_raw']));
data = data.data;

% load recon
elec = load(fullfile(datdir_recon, [sid '_elec_mni_frv']));
elec = elec.elec_mni_frv;

% match channel labels in data to channel labels in recon
[~, idx_data, idx_elec] = intersect(lower(data.label), lower(elec.label));

cfg = [];
cfg.channel = idx_data; % subselect data channels

data = ft_selectdata(cfg, data);

disp(' ');
disp(['channels in recon: ' num2str(length(elec.label))]);
disp(['matching channels in data: ' num2str(length(data.label))]);

% add elec info to data structure and match labels exactly
data.elec.label = elec.label(idx_elec); 
data.elec.elecpos = elec.elecpos(idx_elec,:);
data.elec.chanpos = elec.chanpos(idx_elec,:);
data.label = data.elec.label;

% save data
disp(' '); disp('Saving data structure...');

save(fullfile(savdir, [sid '_data_sync']), 'data', '-v7.3');

clearvars -except data

end
