% Convenience wrapper spanning pipeline 7 and pipeline 8 for one completed run.
%
% Recommended separated entry points:
%   scripts/script_run_one_cfg_bold_reskoopnet_postprocessing.m
%   scripts/script_run_one_bold_cross_modal_coupling.m

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

%% -------------------- user settings --------------------
result_dir = ['E:\autodl_results_local\bold_wsl\e10gb1\mlp\outputs\', ...
    'mlp_obs_bold_wsl_20260423_e10gb1_projected_kv_HP_svd100'];

dataset_stem = 'e10gb1';
dataset_id = '';
observable_mode = 'HP_svd100';

density_sources = struct([]);
density_sources(1).name = 'blp_event_density';
density_sources(1).type = 'event_density';
density_sources(1).file = fullfile( ...
    io_project.get_pipeline_stage_dir(io_project.get_project_processed_root(), dataset_stem, 2, 'event_density'), ...
    sprintf('%s_event_density_2s.mat', dataset_stem));

max_lag_sec = 10;
border_mask_sec = 10;
top_n = 5;

activation_value_mode = 'abs';    % 'abs' or 'real'
activation_slice_list = 1:20;
activation_tiles_per_row = 10;
activation_highlight_core_rois = false;
activation_annotate_regions = false;
activation_annotation_exact_names = {'ParabrachialN', 'Brainstem', 'HP', 'Tha', 'LGN'};
activation_annotation_contains_tokens = {'V1'};
activation_annotation_repeat_each_slice = false;
activation_outline_regions = true;
activation_outline_exact_names = {'ParabrachialN', 'Brainstem', 'HP', 'Tha', 'LGN'};
activation_outline_contains_tokens = {'V1'};

%% -------------------- run combined pipeline 7 + 8 --------------------
warning(['script_run_one_bold_reskoopnet_post_xcorr_activation_maps is now a ', ...
    'convenience wrapper spanning pipeline 7 and pipeline 8. Use the ', ...
    'separated entry points for the canonical workflow.']);

if exist(result_dir, 'dir') ~= 7
    error('result_dir does not exist: %s', result_dir);
end

[~, run_name] = fileparts(result_dir);
processed_root = io_project.get_project_processed_root();
if isempty(dataset_id)
    dataset_id = io_project.get_dataset_id_from_stem(dataset_stem);
end

run_info = struct();
run_info.dataset_stem = dataset_stem;
run_info.dataset_id = dataset_id;
run_info.run_name = run_name;
run_info.output_dir = result_dir;
run_info.observable_mode = observable_mode;

post_params = build_bold_reskoopnet_postprocessing_params();
post_params.processed_root = processed_root;
post_params.skip_existing = false;
post_params.headless = true;
post_params.close_figures = true;
post_result = run_one_bold_reskoopnet_postprocessing_core(run_info, post_params);

p8_params = build_bold_cross_modal_coupling_params();
p8_params.processed_root = processed_root;
p8_params.density_sources = density_sources;
p8_params.skip_existing_xcorr = false;
p8_params.load_existing_xcorr = true;
p8_params.xcorr.max_lag_sec = max_lag_sec;
p8_params.xcorr.border_mask_sec = border_mask_sec;
p8_params.xcorr.top_n = top_n;
p8_params.activation.top_n = top_n;
p8_params.activation.value_mode = activation_value_mode;
p8_params.activation.slice_list = activation_slice_list;
p8_params.activation.tiles_per_row = activation_tiles_per_row;
p8_params.activation.highlight_core_rois = logical(activation_highlight_core_rois);
p8_params.activation.annotate_regions = logical(activation_annotate_regions);
p8_params.activation.annotation_exact_names = activation_annotation_exact_names;
p8_params.activation.annotation_contains_tokens = activation_annotation_contains_tokens;
p8_params.activation.annotation_repeat_each_slice = logical(activation_annotation_repeat_each_slice);
p8_params.activation.outline_regions = logical(activation_outline_regions);
p8_params.activation.outline_exact_names = activation_outline_exact_names;
p8_params.activation.outline_contains_tokens = activation_outline_contains_tokens;
p8_params.activation.skip_existing = false;
p8_params.roi_summary.top_n = top_n;
p8_params.roi_summary.highlight_exact_names = activation_outline_exact_names;
p8_params.roi_summary.highlight_contains_tokens = activation_outline_contains_tokens;
p8_params.roi_summary.skip_existing = false;

one_result = run_one_bold_cross_modal_coupling_core(post_result.post_file, p8_params);
xcorr_out = one_result.xcorr_out;
act_result = one_result.act_result;

fprintf('\nFinished combined pipeline 7 + 8 convenience run.\n');
fprintf('Post MAT:\n  %s\n', post_result.post_file);
fprintf('XCORR MAT:\n  %s\n', one_result.xcorr_mat_file);
fprintf('Activation maps:\n  %s\n', one_result.activation_dir);
fprintf('ROI summaries:\n  %s\n', one_result.roi_summary_dir);
