function ieeg_4_timestamp_wmdms(sid)
% IEEG_4_TIMESTAMP_WMDMS - timestamp the data to mark onsets of trial 
% start times, end times, and/or other events recorded in the photodiode
% signal at the sampling rate of the data. Read behavioral data and save as
% trialinfo variable in standard FieldTrip format. Interactive. This script
% is specific to the WM-DMS task.
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Inputs:
% sid = subject ID (e.g., 'NM01')
%
% Example:
% ieeg_4_timestamp_wmdms('NM01')
%
% Copyright (c) 2022-2025
% EL Johnson, PhD

clearvars -except sid

% set directories
pth = fullfile(pwd, sid);
datdir = fullfile(pth, 'Preprocessing', 'WM-DMS'); % output of ieeg_3_select_task
savdir = datdir;

% load data
data = load(fullfile(datdir, [sid '_data']));
data = data.data;

% plot photodiode signal
photo = rmfield(data, 'events');
photo.label = {'photo'};
photo.trial{1} = data.events{1}.photo;
clear data

cfg = [];
cfg.viewmode = 'vertical';
cfg.blocksize = 30;

ft_databrowser(cfg, photo);

% clean up photodiode signal
disp(' '); invert = input('Enter 1 if the photodiode signal is inverted (up is down): ');

if invert == 1
    photo.trial{1} = photo.trial{1}*(-1);
end

ph = photo.trial{1};

ph_diff = (abs(max(ph))-abs(min(ph))) / 3;
ph_min = min(ph) + ph_diff;
ph_max = max(ph) - ph_diff;
if ph_min > 0
    ph(ph <= ph_min) = 0;
    ph(ph > 0) = 500;
else
    ph(ph >= ph_max) = 500;
    ph(ph < 500) = 0;   
end

% detect onsets as photodiode signal up-changes
e = 0;
for t = 1:length(ph) - 1
    if ph(t+1) == 500 && ph(t) == 0
        e = e + 1;
        onsets(e) = t + 1;
    end
end
clear ph

% read task data (text file)
fid = dir(fullfile(pth, [sid '_trial*']));

tmp = fopen(fullfile(pth, fid.name), 'r');
dat = fscanf(tmp, '%f %f %f %f %f %f %f %f %f %f %f', [13 Inf]);
fclose(tmp);

dat = dat';
trialinfo = dat(1:11,:)';

% align onsets with task data
star = mod(dat(:,3),2) == 1; % identify star trials
clear dat

c = 3; % first 3 up-changes mark task start
idx = zeros(size(onsets));
for x = 1:length(star)
    if star(x) == 0
        c = c + 1;
        idx(c) = 1; 
    elseif star(x) == 1 % 3 onsets/star trial
        c = c + 1;
        idx(c) = 1;
        c = c + 2;
    end
end

onsets = onsets(logical(idx));

% check
disp(' '); disp(['Identified ' num2str(length(onsets)) ' trials...'])
if length(onsets) ~= 128
    error('There should be 128. Mark manually using ft_databrowser.')
end

% set trial start and end times
trl = zeros(length(onsets),3); % initialize
trl(:,1) = onsets; % start times, including 1-s pretrial buffer by design
trl(:,2) = trl(:,1) + (8.5.*photo.fsample); % end times: prestim --> delay 2 with 1-s buffers
trl(:,3) = -1; % start time in s: overwrite with -1 for 1-s pretrial buffer

% check to prevent error if data file ends before final time by making
% final time end of data file
if trl(end,2) > photo.time{1}(end)
   trl(end,2) = photo.time{1}(end);  
end

srate = photo.fsample; % save sampling rate with data for matching later
clear photo

% save
disp('Saving timestamps and trial info...');

save(fullfile(savdir, [sid '_timestamp']), 'trl', 'trialinfo', 'srate');
