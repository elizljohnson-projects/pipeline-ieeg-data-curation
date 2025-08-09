function ieeg_0_read_scans(sid)
% IEEG_0_READ_SCANS - read MRI (pre-op) and CT scans, process MRI using 
% FreeSurfer, fuse MRI and CT scans, and save scans as nifti files. 
% Interactive.
% 
% Analysis described in: Stolk, A, Griffin, S, van der Meij, R et al. 
% Integrated analysis of anatomical and electrophysiological human
% intracranial data. Nature Protocols 13 (2018). 
% https://doi.org/10.1038/s41596-018-0009-6
%
% Guide to fiducials: https://people.cas.sc.edu/rorden/anatomy/na_ac.html
%
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Ensure FreeSurfer home directory is correct in top section of function.
% 
% FreeSurfer steps take hours to days. Use of a high-performance compute
% cluster or similar is recommended.
%
% Inputs:
% sid = subject ID (e.g., 'NM01')
%
% Example:
% ieeg_0_read_scans('NM01')
%
% Copyright (c) 2022-2025
% EL Johnson, PhD

clearvars -except sid

% path to freesurfer (replace with your path to freesurfer home directory)
fshome = '/hpc/software/freesurfer/7.1';

%% set directories

pth = fullfile(pwd, sid);
datdir_scans = fullfile(pth, 'Scans');
savdir = fullfile(pth, 'Recon');
tmp = datetime('today'); % label save folder by date
savdir_recon = fullfile(savdir, ...
    strcat('Recon_', string((month(tmp, 'shortname'))), '_', string(year(tmp))), ...
    'FT_Pipeline');
mkdir(fullfile(savdir_recon, 'Scans'));
mkdir(fullfile(savdir_recon, 'Surfaces'));

% add path to subfunctions
addpath(fullfile(pwd, 'recon_subfunctions'));

clear tmp

%% select scans (1st dicom file in selected folder)
% for more than one MRI and/or CT folder, view each series and keep only
% the selected series folder in the MR/CT folder
%
% search_dicomseries(fullfile(datdir_scans, 'MR'));
% search_dicomseries(fullfile(datdir_scans, 'CT'));
%
% MRI tips:
%   - pick the series with the highest number of dicoms 
%   - stick with series in folders containing 'T1 3D MPRAGE'
%   - avoid series in folders labeled T2, flair, or contrast
%   - avoid contrast scans (fluid appears white)
% CT tips:
%   - pick the series with the highest number of dicoms 
%   - check in DICOM viewer to ensure the image is clear 
%   - ability to see electrodes in white over gray, non-stealth brain bone

dir_mr = dir(fullfile(datdir_scans, 'MR'));
tmp = dir(fullfile(datdir_scans, 'MR', dir_mr(3).name));
fid_mr = fullfile(datdir_scans, 'MR', dir_mr(3).name, tmp(3).name);
clear tmp

dir_ct = dir(fullfile(datdir_scans, 'CT'));
tmp = dir(fullfile(datdir_scans, 'CT', dir_ct(3).name));
fid_ct = fullfile(datdir_scans, 'CT', dir_ct(3).name, tmp(3).name);
clear tmp

%% read and preprocess MRI

% read scan
mri = ft_read_mri(fid_mr);

% determine left/right hemisphere
mri = ft_determine_coordsys(mri);
    % right hemisphere is identified as the hemisphere having larger values 
    % along the left-right axis

% set coordinate system and fiducials
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'acpc';

mri_acpc = ft_volumerealign(cfg, mri);

% save as nifti
disp(' '); disp('Saving preprocessed MRI...');

cfg = [];
cfg.filename = fullfile(savdir_recon, 'Scans', [sid '_MR_acpc']);
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';

ft_volumewrite(cfg, mri_acpc);

mrfile = [cfg.filename '.nii'];

clear mri*

%% process MRI using FreeSurfer (takes hours to days)

disp(' '); disp('Processing MRI using FreeSurfer...');

% setup freesurfer environment
addpath(fullfile(fshome, 'matlab'));
setenv('FREESURFER_HOME', fshome);
PATH = getenv('PATH');
setenv('PATH', [PATH ':' fshome '/bin']);

% initialize nifti file, to be populated and saved in next step
tmp_mrfile = fullfile(datdir_scans, 'MR', [sid '_MR.nii']); 
mkdir(fullfile(savdir, 'freesurfer'));

system(['export FREESURFER_HOME=' fshome '; ' ...
    'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
    'mri_convert -c -oc 0 0 0 ' mrfile ' ' tmp_mrfile '; ' ...
    'recon-all -i ' tmp_mrfile ' -s ' 'freesurfer' ' -sd ' ...
    fullfile(savdir, 'freesurfer') ' -all']);
    % if you get 'flag' error, append with '-cw256'

clear tmp*

% load processed MRI
fsmri_acpc = ft_read_mri(fullfile(savdir, 'freesurfer', 'freesurfer', ...
    'mri', 'T1.mgz'));

% save as nifti
disp(' '); disp('Saving post-processed MRI...');

cfg = [];
cfg.filename = fullfile(savdir_recon, 'Scans', [sid '_fsMR_acpc']); 
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';

ft_volumewrite(cfg, fsmri_acpc);

%% generate surface meshes and hulls from processed MRI (takes hours)

disp(' '); disp('Preparing surface meshes and hulls...');

% prepare hulls 
cfg = [];
cfg.method = 'cortexhull';
cfg.fshome = fshome;

cfg.headshape = fullfile(datdir_scans, 'lh.pial.T1');
hull_lh = ft_prepare_mesh(cfg);

cfg.headshape = fullfile(datdir_scans, 'rh.pial.T1');
hull_rh = ft_prepare_mesh(cfg);

% prepare meshes
cortex_lh = ft_read_headshape(fullfile(datdir_scans, 'lh.pial.T1'));
cortex_rh = ft_read_headshape(fullfile(datdir_scans, 'rh.pial.T1'));

% save
disp(' '); disp('Saving surface meshes and hulls...');

save(fullfile(savdir_recon, 'Surfaces', [sid '_hull_lh']), 'hull_lh');
save(fullfile(savdir_recon, 'Surfaces', [sid '_hull_rh']), 'hull_rh');
save(fullfile(savdir_recon, 'Surfaces', [sid '_cortex_lh']), 'cortex_lh');
save(fullfile(savdir_recon, 'Surfaces', [sid '_cortex_rh']), 'cortex_rh');

%% read and process CT scan

% read scan
ct = ft_read_mri(fid_ct);

% determine left/right hemisphere
ct = ft_determine_coordsys(ct);
    % right hemisphere is identified as the hemisphere having larger values 
    % along the left-right axis
    
% set coordinate system and fiducials
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'acpc';

ct_acpc = ft_volumerealign(cfg, ct);

% save as nifti
disp(' '); disp('Saving preprocessed CT scan...');

cfg = [];
cfg.filename = fullfile(savdir_recon, 'Scans', [sid '_CT_acpc']);
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';

ft_volumewrite(cfg, ct_acpc);

%% fuse MRI and CT

cfg = [];
cfg.method = 'spm';
cfg.spmversion = 'spm12';
cfg.coordsys = 'acpc';
cfg.viewresult = 'yes';

ct_acpc_f = ft_volumerealign(cfg, ct_acpc, fsmri_acpc);
ct_acpc_f.coordsys = 'acpc';

clear ct_acpc

% save figure of result for reference
print(fullfile(savdir_recon, 'Scans', [sid '_CT_fsMR_fusion.png']), '-dpng');
  
% save realigned CT as nifti
disp(' '); disp('Saving post-processed CT scan...');

cfg = [];
cfg.filename = fullfile(savdir_recon, 'Scans', [sid '_CT_acpc_f']);
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';

ft_volumewrite(cfg, ct_acpc_f);

close all

end
