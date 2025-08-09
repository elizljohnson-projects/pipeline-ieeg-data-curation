function ieeg_6_reref_data(sid, task, ch_fid)
% IEEG_6_REREF_DATA - re-reference clean data using bipolar montage. The
% channel pair lists are input by the user, and should follow standard
% conventions (e.g., sEEG deep-surface, ECoG anterior-posterior). Channel 
% pairs should be adjacent (e.g., LOF1, LOF2), follow anatomical boundaries
% (e.g., same major gyrus), and exclude white matter-white matter pairs.
% The list should be saved as a .csv file in the ~/Preprocessing folder.
% 
% Identify and remove and residual trials and re-referenced channels 
% containing artifacts. Interactive.
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Inputs:
% sid = subject ID (e.g., 'NM01')
% task = task name (e.g., 'WM-DMS')
% ch_fid = filename of channel pair list with 2 columns: channels to 
%   subtract from in column 1, channels to subtract in column 2
%
% Example:
% ieeg_6_reref_data('NM01', 'WM-DMS', 'NM01_ch_reref.csv')
%
% Copyright (c) 2022
% EL Johnson, PhD

clearvars -except sid task ch_fid

% set directories
pth = fullfile(pwd, sid);
datdir = fullfile(pth, 'Preprocessing', task); % output of ieeg_5_clean_data
datdir_montage = fullfile(pth, 'Preprocessing');
savdir_montage = datdir_montage;
savdir_data = fullfile(pth, task);
mkdir(savdir_data);

% load data
data = load(fullfile(datdir, [sid '_data_clean']));
data = data.data;

% load channel list
ch_list = readtable(fullfile(datdir_montage, ch_fid));

% find channels in data iteratively in case of channel removal
ch1 = cell(size(ch_list,1), 1);
ch2 = ch1;
idx1 = nan(size(ch1));
idx2 = idx1;
for c = 1:length(ch1)
    [ch1{c}, idx1(c)] = intersect(data.label, ch_list{c,1});
    if ~empty(ch1{c})
        [ch2{c}, idx2{c}] = intersect(data.label, ch_list{c,2});
    end
end
ch1(empty(ch1)) = [];
ch2(empty(ch2)) = [];
idx1(empty(idx1)) = [];
idx2(empty(idx2)) = [];

% save montage
save(fullfile(savdir_montage, [sid '_bpm_reref']), 'ch1', 'ch2', 'idx*');

% re-reference
cfg = [];
cfg.channel = ch1;

data = ft_selectdata(cfg, data);

for r = 1:length(data.trial)
    data.trial{r} = data.trial{r}(idx1,:) - data.trial{r}(idx2,:);
end

% overwrite channel labels and compute virtual x-y-x coordinates as 
% channel pair means

data.label = cell(size(data.label));
for c = 1:length(ch1)
    data.label{c} = strcat(ch1{c}, '-', ch2{c}); % new label
end

% virtual coordinates
data.elec.label = data.label;
data.elec.elecpos = (data.elec.elecpos(idx1,:) + data.elec.elecpos(idx2,:)) ./ 2;
data.elec.chanpos = (data.elec.chanpos(idx1,:) + data.elec.chanpos(idx2,:)) ./ 2;

% identify any remaining artifacts
cfg = [];
cfg.viewmode = 'vertical';
cfg.allowoverlap = 'yes';

disp(' '); disp(' ');
disp('Mark any remaining noisy epochs for removal, then close plotting window.');
disp('Also write down any remaining noisy channels to be removed.')

artfct = ft_databrowser(cfg, data);

% remove marked trials
cfg = [];
cfg.artfctdef.visual.artifact = artfct.artfctdef.visual.artifact;
cfg.allowoverlap = 'yes';

data = ft_rejectartifact(cfg, data);
clear artfct

% remove any remaining bad channels
disp(' '); bad_chans = input('Enter channels to remove in {}: ');

cfg = [];
cfg.channel = setxor(data.label, bad_chans);

data = ft_selectdata(cfg, data);

% save
disp(' '); disp('Saving data, ready to analyze...');

save(fullfile(savdir, sid), 'data');

end
