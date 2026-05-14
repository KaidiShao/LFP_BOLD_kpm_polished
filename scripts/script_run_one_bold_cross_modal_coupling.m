% Convenience wrapper for running pipeline 8 on one existing BOLD_POST file.
%
% Recommended canonical entry:
%   scripts/script_run_one_cfg_bold_cross_modal_coupling.m

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

%% -------------------- user settings --------------------
if ~exist('bold_post_file', 'var') || isempty(bold_post_file)
    error(['Set bold_post_file before running this script. Example:\n', ...
        'bold_post_file = ''E:\\DataPons_processed\\e10gb1\\...\\mat\\..._bold_post.mat'';']);
end

params = build_bold_cross_modal_coupling_params();
params.headless = true;
params.close_figures = true;
params.load_existing_xcorr = true;
params.skip_existing_xcorr = false;
params.activation.skip_existing = false;

if exist('density_sources', 'var') && ~isempty(density_sources)
    params.density_sources = density_sources;
end
if exist('force_recompute', 'var') && ~isempty(force_recompute)
    params.skip_existing_xcorr = ~logical(force_recompute);
    params.activation.skip_existing = ~logical(force_recompute);
end
if exist('max_lag_sec', 'var') && ~isempty(max_lag_sec)
    params.xcorr.max_lag_sec = max_lag_sec;
end
if exist('border_mask_sec', 'var') && ~isempty(border_mask_sec)
    params.xcorr.border_mask_sec = border_mask_sec;
end
if exist('top_n', 'var') && ~isempty(top_n)
    params.xcorr.top_n = top_n;
    params.activation.top_n = top_n;
end
if exist('activation_value_mode', 'var') && ~isempty(activation_value_mode)
    params.activation.value_mode = char(string(activation_value_mode));
end
if exist('activation_slice_list', 'var') && ~isempty(activation_slice_list)
    params.activation.slice_list = activation_slice_list;
end
if exist('activation_tiles_per_row', 'var') && ~isempty(activation_tiles_per_row)
    params.activation.tiles_per_row = activation_tiles_per_row;
end
if exist('activation_feature_reduce', 'var') && ~isempty(activation_feature_reduce)
    params.activation.feature_reduce = char(string(activation_feature_reduce));
end

one_result = run_one_bold_cross_modal_coupling_core(bold_post_file, params);
xcorr_out = one_result.xcorr_out;
act_result = one_result.act_result;

fprintf('\nFinished pipeline 8 for one BOLD ResKoopNet run.\n');
fprintf('XCORR status: %s\n', one_result.xcorr_status);
fprintf('Activation status: %s\n', one_result.activation_status);
fprintf('XCORR MAT:\n  %s\n', one_result.xcorr_mat_file);
fprintf('Activation maps:\n  %s\n', one_result.activation_dir);
