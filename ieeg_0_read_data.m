function ieeg_0_read_data(sid, fid)
% IEEG_0_READ_DATA - read raw iEEG data and save as standardized MATLAB 
% structure in FieldTrip format.
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Inputs:
% sid = subject ID (e.g., 'NM01')
% fid = raw iEEG data filename (e.g., 'NM01.edf') or folder name (e.g.,
%   '2021-12-02_14-59-09')
%
% Example:
% ieeg_0_read_data('NM01', '2021-12-02_14-59-09')
%
% Copyright (c) 2022-2025
% EL Johnson, PhD

clearvars -except sid fid

% set directories
pth = fullfile(pwd, sid);
savdir = fullfile(pth, 'Preprocessing');
mkdir(savdir)

% initialize data structure
data = [];

% read data and put in structure
disp(' '); disp('Reading data file or folder...'); disp(' ');

% if file, read file (e.g., Nihon Kohden/Natus data)
if contains(fid, '.')
    datdir = pth;

    % read header and data
    tmp_hdr = ft_read_header(fullfile(datdir, fid));
    tmp_data = ft_read_data(fullfile(datdir, fid));

    % populate structure
    data.fsample = tmp_hdr.Fs;
    data.label = tmp_hdr.label;
    data.trial{1} = tmp_data;
    data.time{1} = linspace(0, length(tmp_data), length(tmp_data));
    clear tmp*

    % clean up channel names
    for c = 1:length(data.label)
        if ~isempty(strfind(data.label{c}, ' ')) % often prepended 'EEG '
            data.label{c} = data.label{c}(strfind(data.label{c}, ' ') + 1:end);
        end
        if ~isempty(strfind(data.label{c}, '-')) % often appended '-Ref'
            data.label{c} = data.label{c}(1:strfind(data.label{c}, '-') - 1);
        end
    end

    % add event info from recording notes
    disp(' ');
    ph = input('Enter photodiode channel (e.g., DC1, DC03): ', 's');
    [~, ph_idx] = intersect(data.label, {ph});
    data.events.photo = data.trial(ph_idx,:);

    % mic = input('Enter microphone channel (e.g., DC2, DC04): ', 's');
    % [~, mic_idx] = intersect(data.label, {mic});
    % data.events.mic = data.trial(mic_idx,:);

% if folder, read by channel (e.g., Neuralynx data)
else
    datdir = fullfile(pth, fid);

    % finish initializing data structure based on sample channel
    if f == 1
        % get list of channel filenames
        ch = dir(fullfile(datdir, 'CSC*')); % for Neuralynx; update for
            % per-channel data with another labeling convention

        % read header
        tmp_hdr = ft_read_header(fullfile(datdir, ch(1).name));

        data.fsample = tmp_hdr.Fs;
        data.label = cell(length(ch), 1);
    end

    % read data
    tmp_data = ft_read_data(fullfile(datdir, ch(1).name));

    data.trial = zeros(length(ch), length(tmp_data));
    data.time = linspace(0, length(tmp_data), length(tmp_data));
    clear tmp*

    % add event info
    ev = dir(fullfile(datdir, 'Events*'));
    data.events = ft_read_event(fullfile(datdir, ev.name));

    % populate structure
    for c = 1:length(ch)
        data.trial(c,:) = ft_read_data(fullfile(datdir, ch(c).name));

        % clean up channel names
        if isempty(strfind(ch(c).name, '_')) % often appended '_000X.ncs'
            data.label{c} = ch(c).name(1:strfind(ch(c).name, '.') - 1);
        else
            data.label{c} = ch(c).name(1:strfind(ch(c).name, '_') - 1);
        end
        if length(data.label{c}) == 4
            data.label{c} = strcat(data.label{c}(1:3), '00', data.label{c}(4));
        elseif length(data.label{c}) == 5
            data.label{c} = strcat(data.label{c}(1:3), '0', data.label{c}(4:5));
        end
    end

    % match channel order to channel numbers (e.g., 1 is first, not 10)
    if str2double(data.label{1}(4:6)) > 1
        [data.label, ch_order] = sort(data.label);
        data.trial = data.trial(ch_order,:);
    end
end

% save channel list for the electrode reconstruction
chan = data.label;
save(fullfile(savdir, [sid '_chan_orig']), 'chan');

% save data
disp(' '); disp('Saving data structure...');

save(fullfile(savdir, [sid '_data_raw']), 'data', '-v7.3');

end
