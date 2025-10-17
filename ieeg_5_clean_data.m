function ieeg_5_clean_data(sid, task)
% IEEG_5_CLEAN_DATA - filter (bandpass 0.1-300 Hz, notch line noise), 
% down-sample (if >1 kHz), and segment continuous data. Identify and remove 
% trials and channels containing seizure activity and artifacts. 
% Interactive.
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Filtering steps take minutes to hours. High-performance computing is
% recommended.
%
% Inputs:
% sid = subject ID (e.g., 'NM01')
% task = task name (e.g., 'WM-DMS')
%
% Example:
% ieeg_5_clean_data('NM01', 'WM-DMS')
%
% Copyright (c) 2022-2025
% EL Johnson, PhD

clearvars -except sid task

% set directories
pth = fullfile(pwd, sid);
datdir = fullfile(pth, 'Preprocessing', task); % output of ieeg_4_timestamp
savdir = datdir;

% load data
data = load(fullfile(datdir, [sid '_data']));
data = data.data;

% load timestamps and trialinfo
tmp = load(fullfile(datdir, [sid '_timestamp']));
trl = tmp.trl;
trialinfo = tmp.trialinfo;
srate = tmp.srate;
clear tmp

% bandpass filter 0.1-300 Hz (takes minutes to hours)
data.trial{1} = eegfilt(data.trial{1}, data.fsample, 0.1, [], 0, 0, 0, 'fir1');
data.trial{1} = eegfilt(data.trial{1}, data.fsample, [], 300, 0, 0, 0, 'fir1');

% notch filter to remove line noise
cfg = [];
cfg.dftfilter = 'yes';
cfg.dftfreq = [60 120 180 240]; % 60 Hz for USA, 50 Hz for Europe
cfg.demean = 'yes';

data = ft_preprocessing(cfg, data);

% down-sample if recorded > 1 kHz
if data.fsample > 1000
    trl = round(trl .* (1000/srate)); % resample timestamps 
    data.trial{1} = double(data.trial{1}); 
    [p,q] = rat(1000/data.fsample);
        
    for e = 1:size(data.trial{1},1)
        tmp = resample(data.trial{1}(e,:), p, q);        
        if e == 1
            tmp_mtx = nan(size(data.trial{1},1), length(tmp));
        end        
        tmp_mtx(e,:) = tmp; 
        clear tmp
    end
    
    data.fsample = 1000;
    data.trial{1} = tmp_mtx;
    data.time{1} = linspace(0, size(tmp_mtx,2), size(tmp_mtx,2));
    clear tmp_mtx    
end

% save filtered data
disp(' '); disp('Saving filtered data...');

save(fullfile(savdir, [sid '_data_filt']), 'data');

% segment continuous data into trials
cfg = [];
cfg.trl = trl;

data = ft_redefinetrial(cfg, data);

data = rmfield(data, 'events');
data.trialinfo = trialinfo; % add trial info
clear tmp trialinfo

% save trial-segmented data
disp(' '); disp('Saving trial-segmented data...');

save(fullfile(savdir, [sid '_data_trl']), 'data');

% identify artifacts
cfg = [];
cfg.viewmode = 'vertical';
cfg.allowoverlap = 'yes';

disp(' '); disp(' ');
disp('Mark noisy epochs for removal, then close plotting window.');
disp('Note: the whole trial will be removed if marked.')
disp('Also write down noisy channels to be removed.')

artfct = ft_databrowser(cfg, data);

% save marked artifacts
save(fullfile(savdir, [sid '_artfct']), 'artfct');

% remove marked trials
cfg = [];
cfg.artfctdef.visual.artifact = artfct.artfctdef.visual.artifact;
cfg.allowoverlap = 'yes';

data = ft_rejectartifact(cfg, data);
clear artfct

% remove bad channels
disp(' '); bad_chans = input('Enter channels to remove in {}: ');

cfg = [];
cfg.channel = setxor(data.label, bad_chans);

data = ft_selectdata(cfg, data);

% save bad channels
save(fullfile(savdir, [sid '_bad_chans']), 'bad_chans');
clear bad_chans

% save clean data
disp(' '); disp('Saving clean data...')

save(fullfile(savdir, [sid '_data_clean']), 'data');

end

