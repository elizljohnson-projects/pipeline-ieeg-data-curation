function ieeg_3_select_task(sid)
% IEEG_3_SELECT_TASK - cut continuous data into task segments and save in
% separate data structures. Interactive.
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Inputs:
% sid = subject ID (e.g., 'NM01')
%
% Example:
% ieeg_3_select_task('NM01')
%
% Copyright (c) 2022-2025
% EL Johnson, PhD

clearvars -except sid

% set directories
pth = fullfile(pwd, sid);
datdir = fullfile(pth, 'Preprocessing'); % output of ieeg_2_sync_data
savdir = datdir;

% load data
data = load(fullfile(datdir, [sid '_data_sync']));
data = data.data;

% segment based on event data

% create event structure for photodiode (and microphone etc.)
events = rmfield(data, 'events');
fn = fieldnames(data.events{1});
events.label = cell(length(fn), 1);
events.trial{1} = nan(length(fn), length(events.time{1}));
for c = 1:length(fn)
    events.label(c) = fn(c);
    events.trial{1}(c,:) = data.events{1}.(fn{c});
end

% plot and manually mark task start and end times
cfg = [];
cfg.viewmode = 'vertical';
cfg.blocksize = 30;

disp(' '); disp(' ');
disp('Mark task start and end times, including some buffer, then close plotting window.');
disp('Note: the photodiode signal may be inverted (up is down).')

times = ft_databrowser(cfg, events);

% segment the data into tasks
task_count = size(times.artfctdef.visual.artifact, 1) / 2;

disp(' '); disp(['Splitting the data into ' num2str(task_count) ' tasks...']);

tasks = cell(task_count, 1);
for t = 1:2:size(times.artfctdef.visual.artifact, 1)
    cfg = [];
    cfg.begsample = times.artfctdef.visual.artifact(t,1);
    cfg.endsample = times.artfctdef.visual.artifact(t+1,1);

    tasks{t} = ft_redefinetrial(cfg, data);
end
clear data times

% save in folder by task name
disp(' '); disp('Saving the data by task name...'); disp(' ');

for t = 1:task_count
    data = tasks{t};
    tmp_name = input(['Enter name of task #' num2str(t) ...
        ' (e.g., WM-DMS, MemDev, Rest): '], 's');
    mkdir(fullfile(savdir, tmp_name))
    save(fullfile(savdir, tmp_name, [sid '_data']), 'data', '-v7.3');
    clear data tmp*
end

end
