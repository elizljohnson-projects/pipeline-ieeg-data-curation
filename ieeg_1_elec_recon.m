function ieeg_1_elec_recon(sid)
% IEEG_1_ELEC_RECON - reconstruct electrode positions in native and MNI 
% sapce and save as standardized MATLAB structure in FieldTrip format. Also
% save Excel table of region probability labels in MNI space, and figures 
% of reconstructions in native and MNI space. Interactive.
%
% Analysis described in: Stolk, A, Griffin, S, van der Meij, R et al.
% Integrated analysis of anatomical and electrophysiological human
% intracranial data. Nature Protocols 13 (2018).
% https://doi.org/10.1038/s41596-018-0009-6
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Download and unzip FieldTrip MNI templates before running the first time: 
% https://drive.google.com/file/d/1qWP-v-ytWxXUuKWydSksg1X1hvZjBi8P/view?usp=sharing
%
% Inputs:
% sid = subject ID (e.g., 'NM01')
%
% Example:
% ieeg_1_elec_recon('NM01')
%
% Copyright (c) 2022-2025
% EL Johnson, PhD

clearvars -except sid

% set directories
pth = fullfile(pwd, sid);
datdir_data = fullfile(pth, 'Preprocessing'); % output of ieeg_0_read_data
tmp = dir(fullfile(pth, 'Recon', 'Recon_*'));
datdir_recon = fullfile(pth, 'Recon', tmp.name, 'FT_Pipeline');
datdir_scans = fullfile(datdir_recon, 'Scans');
    % ouput of ieeg_0_read_scans
savdir_scans = datdir_scans;
savdir_elec = fullfile(datdir_recon, 'Electrodes');
mkdir(savdir_elec);
clear tmp

savdir_native = fullfile(pth, 'Recon', 'Native_Brain_2D_Recons');
savdir_mni = fullfile(pth, 'Recon', 'Standardized_Brain_2D_Recons');
mkdir(savdir_native);
mkdir(savdir_mni);

xldir = fullfile(pwd, 'recon_subfunctions', 'xlwrite'); % Excel functions

% add paths to subfunctions and templates
addpath(fullfile(pwd, 'recon_subfunctions'));
addpath(fullfile(pwd, 'ft_templates'));

% load processed scans and set coordinate system
fsmri_acpc = ft_read_mri(fullfile(datdir_scans, [sid '_fsMR_acpc.nii']));
fsmri_acpc.coordsys = 'acpc';

ct_acpc_f = ft_read_mri(fullfile(datdir_scans, [sid '_CT_acpc_f.nii']));
ct_acpc_f.coordsys = 'acpc';

%% manually place electrodes

% load channel list from raw iEEG data
chan = load(fullfile(datdir_data, [sid '_chan_orig']));

cfg = [];
cfg.channel = chan.chan;
clear chan

% 1. Orthoplot viewing options:
%    a. use the left mouse button to navigate the image, or
%    b. use the arrow keys to increase or decrease the slice number by one
% 2. Orthoplot placement options:
%    a. click an electrode label in the list to assign it to the crosshair
%       location, or
%    b. doubleclick a previously assigned electrode label to remove its
%       marker
% 3. To finalize, close the window or press q on the keyboard

elec_acpc_f = ft_electrodeplacement(cfg, ct_acpc_f, fsmri_acpc);

% Tips:
%   ECoG grid numbering: check patient notes
%   sEEG depth numbering: deepest electrode is always #1

clear ct_acpc_f

%% describe electrode types

cfg = [];
cfg.detchantype = 'auto';

elec_acpc_f = describe_elec(cfg, elec_acpc_f);

disp(' '); disp('Review electrode description:'); disp(elec_acpc_f);
disp(' '); elec_check = input('Enter 1 to manually adjust the description: ');

% correct electrode types if needed
if elec_check == 1
    disp(' '); ecog = input('Enter 1 if the dataset contains ECoG grids/strips: ');
    seeg = input('Enter 1 if the dataset contains sEEG depth tracks: ');

    elec_acpc_f.surface = {};
    elec_acpc_f.depths = {};
    
    disp(' ');
    if ecog == 1
        grid_hold = input('Enter number of separate grids + strips (e.g., 2): ');
        grids = cell(grid_hold,1);
        grid = [];
        for x = 1:grid_hold
            tmp = input(['Enter channel numbers of grid/strip #' num2str(x) ' (e.g., 1:32): ']);
            grid = cat(2, grid, tmp);
            grids{x} = elec_acpc_f.label(tmp);
            clear tmp
        end
        for e = grid
            elec_acpc_f.chantype{e} = 'ieeg_grid';
            % strips and grids treated the same way, so use the same label for
            % manual overwrite
        end
        elec_acpc_f.surface{1} = elec_acpc_f.label(grid);
    end
    
    if seeg == 1
        depth_hold = input('Enter number of separate depth tracks (e.g., 2): ');
        depths = cell(depth_hold,1);
        depth = [];
        for x = 1:depth_hold
            tmp = input(['Enter channel numbers of depth track #' num2str(x) ' (e.g., 39:48): ']);
            depth = cat(2, depth, tmp);
            depths{x} = elec_acpc_f.label(tmp);
            clear tmp
        end
        for e = depth
            elec_acpc_f.chantype{e} = 'ieeg_depth';
        end
        elec_acpc_f.depths{1} = elec_acpc_f.label(depth);
    end    
else
    % pull electrode info from description
    if ~isempty(elec_acpc_f.surface)
        ecog = 1;
        grids = elec_acpc_f.surface;
    else
        ecog = 0;
    end
    if ~isempty(elec_acpc_f.depths)
        seeg = 1;
    else
        seeg = 0;
    end
end

% save electrode file
disp(' '); disp('Saving raw electrode file...');

save(fullfile(savdir_elec, [sid '_elec_acpc_f']), 'elec_acpc_f');

%% ECoG grid/strip projection via Dykstra method

% load processed MRI surface meshes
tmp = load(fullfile(datdir_recon, 'Surfaces', [sid '_cortex_lh']));
cortex_lh = tmp.cortex_lh;
clear tmp

tmp = load(fullfile(datdir_recon, 'Surfaces', [sid '_cortex_rh']));
cortex_rh = tmp.cortex_rh;
clear tmp

if ecog == 1
    % load processed MRI hulls
    tmp = load(fullfile(datdir_recon, 'Surfaces', [sid '_hull_lh']));
    hull_lh = tmp.hull_lh;
    clear tmp

    tmp = load(fullfile(datdir_recon, 'Surfaces', [sid '_hull_rh']));
    hull_rh = tmp.hull_rh;
    clear tmp

    ecog_ch = cell(length(grids, 1));
    for x = 1:length(grids)
        cfg = [];
        cfg.channel = grids{x};
        cfg.elec = elec_acpc_f;
        cfg.casesensitive = 'yes';
        cfg.method = 'headshape';
        if strcmp(cfg.elec.chanside{cfg.channel(1)}, 'left')
            cfg.headshape = hull_rh;
        elseif strcmp(cfg.elec.chanside{cfg.channel(1)}, 'right')
            cfg.headshape = hull_lh;
        end
        cfg.warp = 'dykstra2012';
        cfg.pairmethod = 'label';
        cfg.isodistance = 'yes';
        cfg.deformweight = 5; % >1 = more cost to deforming the grid =
            % more likely to maintain its structure

        ecog_ch{x} = ft_electroderealign(cfg);
        ecog_ch{x}.unit = 'mm';
    end

    % create new structure for projected ECoG channels in the chanpos
    % field, with no change to the elecpos field (or sEEG channels if any)
    elec_acpc_fr = elec_acpc_f;
    for x = 1:length(grids)
        elec_acpc_fr.chanpos(grids{x},:) = ecog_ch{x}.chanpos;
    end

    % view projection and adjust if needed
    figure;
    if any(strcmp(unique(elec_acpc_f.chanside(match_str(elec_acpc_f.chantype, ...
            {'ieeg_strip', 'ieeg_grid'}))), 'right')) && ...
            any(strcmp(unique(elec_acpc_f.chanside(match_str(elec_acpc_f.chantype, ...
            {'ieeg_strip', 'ieeg_grid'}))), 'left'))
        ft_plot_mesh(cortex_lh, 'facecolor', [0.781 0.762 0.664], ...
            'EdgeColor', 'none', 'label', 'label');
        ft_plot_mesh(cortex_rh, 'facecolor', [0.781 0.762 0.664], ...
            'EdgeColor', 'none', 'label', 'label');
        view([-55 10]);
    elseif strcmp(unique(elec_acpc_f.chanside(match_str(elec_acpc_f.chantype, ...
            {'ieeg_strip', 'ieeg_grid'}))), 'left')
        ft_plot_mesh(cortex_lh, 'facecolor', [0.781 0.762 0.664], ...
            'EdgeColor', 'none', 'label', 'label');
        view([-115 0]);
    elseif strcmp(unique(elec_acpc_f.chanside(match_str(elec_acpc_f.chantype, ...
            {'ieeg_strip', 'ieeg_grid'}))), 'right')
        ft_plot_mesh(cortex_rh, 'facecolor', [0.781 0.762 0.664], ...
            'EdgeColor', 'none', 'label', 'label');
        view([90 0]);
    end
    lighting gouraud; camlight; material dull;
    ft_plot_sens(elec_acpc_fr ,'elecshape', 'sphere', 'label', 'label');

    % save electrode file
    disp(' '); disp('Saving electrode file with ECoG channels projected...');

    save(fullfile(savdir_elec, [sid '_elec_acpc_fr']), 'elec_acpc_fr');

    clear elec_acpc_f ecog_*
end

%% convert native coordinates to MNI using volume based normalization

disp(' '); disp('Converting native coordinates to MNI coordinates...');

cfg = [];
cfg.nonlinear = 'yes';
cfg.template = fullfile(pwd, 'ft_templates', 'single_subj_T1_1mm.nii'); % MNI template
cfg.spmversion = 'spm12';
cfg.spmmethod = 'new';

fsmri_mni_v = ft_volumenormalise(cfg, fsmri_acpc);

% warp electrode positions to MNI template
if exist('elec_acpc_fr', 'var')
    elec_nat = elec_acpc_fr;
else
    elec_nat = elec_acpc_f;
end
elec_mni_frv = elec_nat;
elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni_v.params, elec_nat.elecpos, ...
    'individual2sn');
elec_mni_frv.chanpos = elec_mni_frv.elecpos;
elec_mni_frv.coordsys = 'mni';

clear elec_acpc*

%% snap ECoG grid/strip channels to MNI template using Dykstra projection

% load MNI surfaces
mni_brain = load(fullfile(pwd, 'ft_templates', 'template_hull'));
surface_mni_lh = load(fullfile(pwd, 'ft_templates', 'surface_pial_left'));
surface_mni_rh = load(fullfile(pwd, 'ft_templates', 'surface_pial_right'));

if ecog == 1   
    ecog_mni = cell(length(grids, 1));
    for x = 1:length(grids)
        cfg = [];
        cfg.channel = grids{x};
        cfg.elec = elec_mni_frv;
        cfg.casesensitive = 'yes';
        cfg.method = 'headshape';
        if strcmp(cfg.elec.chanside{cfg.channel(1)}, 'left')
            cfg.headshape = mni_brain.hull_lh;
        elseif strcmp(cfg.elec.chanside{cfg.channel(1)}, 'right')
            cfg.headshape = mni_brain.hull_rh;
        end
        cfg.warp = 'dykstra2012';
        cfg.pairmethod = 'label';
        cfg.deformweight = 3 ; % >1 = more cost to deforming the grid =
            % more likely to maintain its structure
        cfg.isodistance = 'yes';

        ecog_mni{x} = ft_electroderealign(cfg);
        ecog_mni{x}.unit = 'mm';
    end

    % overwrite projected ECoG channels in the chanpos field, with no 
    % change to the elecpos field (or sEEG channels if any)
    for x = 1:length(grids)
        elec_mni_frv.chanpos(grids{x},:) = ecog_mni{x}.chanpos;
    end

    % view projection
    figure;
    if any(strcmp(unique(elec_mni_frv.chanside(match_str(elec_mni_frv.chantype, ...
            {'ieeg_strip', 'ieeg_grid'}))), 'right')) && ...
            any(strcmp(unique(elec_mni_frv.chanside(match_str(elec_mni_frv.chantype, ...
            {'ieeg_strip', 'ieeg_grid'}))), 'left'))
        ft_plot_mesh(surface_mni_lh.mesh, 'facecolor', [0.781 0.762 0.664], ...
            'EdgeColor', 'none', 'label', 'label');
        ft_plot_mesh(surface_mni_rh.mesh, 'facecolor', [0.781 0.762 0.664], ...
            'EdgeColor', 'none', 'label', 'label');
        view([-55 10]);
    elseif strcmp(unique(elec_mni_frv.chanside(match_str(elec_mni_frv.chantype, ...
            {'ieeg_strip', 'ieeg_grid'}))), 'left')
        ft_plot_mesh(surface_mni_lh.mesh, 'facecolor', [0.781 0.762 0.664], ...
            'EdgeColor', 'none', 'label', 'label');
        view([-115 0]);
    elseif strcmp(unique(elec_mni_frv.chanside(match_str(elec_mni_frv.chantype, ...
            {'ieeg_strip', 'ieeg_grid'}))), 'right')
        ft_plot_mesh(surface_mni_rh.mesh, 'facecolor', [0.781 0.762 0.664], ...
            'EdgeColor', 'none', 'label', 'label');
        view([90 0]);
    end
    lighting gouraud; camlight; material dull;
    ft_plot_sens(elec_mni_frv ,'elecshape', 'sphere', 'label', 'label');

    clear ecog_*
end

%% save normalied MRI and electrode file

disp(' '); disp('Saving MNI-normalized MRI and electrode file...');

% save MRI
cfg = [];
cfg.filename = fullfile(savdir_scans, [sid  '_fsMR_mni_v']);
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';

ft_volumewrite(cfg, fsmri_mni_v);

% save electrode file
save(fullfile(savdir_elec, [sid '_elec_mni_frv']), 'elec_mni_frv');

%% generate table of regional probability labels

disp(' '); disp('Generating table of regional probability labels...');

generate_electable_v3(fullfile(pth, [sid '_Elec_Notes']), 'elec_mni', ...
    elec_mni_frv, 'xldir', xldir, 'fsdir', fullfile(pth, 'Recon', ...
    'freesurfer', 'freesurfer'), 'elec_nat', elec_nat);

%% generate 2D recons

% native space
disp(' '); disp('Generating recon plots in native space...');

cfg = [];
cfg.ptid = sid;
cfg.cortex_rh = cortex_rh;
cfg.cortex_lh = cortex_lh;
cfg.mri = fsmri_acpc;
cfg.method = 'ind';
cfg.imagesdir = fullfile(savdir_native);

if seeg == 1 && ecog == 0
    cfg.electype = 'depth';
elseif seeg == 0 && ecog == 1
    cfg.electype = 'surface';
else
    cfg.electype = 'all';
end

generate_2d_recons(cfg, elec_nat);

clear cortex*

% MNI space
disp(' '); disp('Generating recon plots in MNI space...');

cfg = [];
cfg.ptid = sid;
cfg.method = 'ind';
cfg.cortex_rh = surface_mni_rh.mesh;
cfg.cortex_lh = surface_mni_lh.mesh;
cfg.hullpos =  [mni_brain.hull_rh.pos; mni_brain.hull_lh.pos];
cfg.mri = ft_read_mri(fullfile(pwd, 'ft_templates', 'single_subj_T1_1mm.nii'));
cfg.imagesdir = fullfile(savdir_mni);

if seeg == 1 && ecog == 0
    cfg.electype = 'depth';
elseif seeg == 0 && ecog == 1
    cfg.electype = 'surface';
else
    cfg.electype = 'all';
end

generate_2d_recons(cfg, elec_mni_frv);

close all

end
