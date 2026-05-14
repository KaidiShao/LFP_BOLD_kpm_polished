% Compute pipeline2-style density outputs for legacy E10gb1 PL event files.
%
% This treats mevt_pl and sevt_pl as parallel legacy event detectors and
% writes their density outputs under:
%   E:\DataPons_processed\e10gb1\pipeline2_legacy_event_density
%
% Optional controls, set before run(...) if needed:
%   legacy_group_name = 'spont';
%   density_bin_sec = 2;
%   density_smooth_sigma_sec = 2;

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

if ~exist('legacy_group_name', 'var') || isempty(legacy_group_name)
    legacy_group_name = 'spont';
end
if ~exist('density_bin_sec', 'var') || isempty(density_bin_sec)
    density_bin_sec = 2;
end
if ~exist('density_smooth_sigma_sec', 'var') || isempty(density_smooth_sigma_sec)
    density_smooth_sigma_sec = density_bin_sec;
end

cfg = cfg_E10gb1();
output_root = io_project.get_project_processed_root();

params = struct();
params.bin_sec = double(density_bin_sec);
params.smooth_sigma_sec = double(density_smooth_sigma_sec);

legacy_specs = struct( ...
    'detector_name', {'mevt_pl', 'sevt_pl'}, ...
    'legacy_variable', {'mevt_pl', 'sevt_pl'});

results = cell(numel(legacy_specs), 1);
for i = 1:numel(legacy_specs)
    params.legacy_variable = legacy_specs(i).legacy_variable;
    source_event_file = fullfile( ...
        cfg.raw_data_root, ...
        char(string(legacy_group_name)), ...
        sprintf('%s_%s_%s.mat', cfg.file_stem, char(string(legacy_group_name)), legacy_specs(i).detector_name));

    results{i} = compute_legacy_blp_event_density( ...
        cfg, ...
        output_root, ...
        source_event_file, ...
        legacy_specs(i).detector_name, ...
        params);

    fprintf('Saved legacy %s event-density result to:\n  %s\n', ...
        legacy_specs(i).detector_name, results{i}.save_file);
end
