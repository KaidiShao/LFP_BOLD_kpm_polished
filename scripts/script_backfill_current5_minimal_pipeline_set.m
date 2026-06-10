% Backfill the current five datasets to the minimum comparison-ready set.
%
% This script is intentionally conservative:
%   - it only runs missing MATLAB-side stages;
%   - it focuses on the minimal cross-dataset comparison set;
%   - it does not train BOLD ResKoopNet models; P4 BOLD training is manual.
%
% Typical dry-run from MATLAB:
%   dry_run = true;
%   backfill_phase = 'all';
%   p3_bold_modes = {'global_svd100','gsvd100_ds','HP_svd100','roi_mean'};
%   minimal_bold_modes = {'global_svd100','gsvd100_ds','HP_svd100','roi_mean'};
%   run('scripts/script_backfill_current5_minimal_pipeline_set.m');
%
% Execute from MATLAB:
%   dry_run = false;
%   backfill_phase = 'all';
%   run('scripts/script_backfill_current5_minimal_pipeline_set.m');
%
% Phases:
%   pre_bold  - P3 BOLD observables/QC
%   blp       - P5 BLP eigenfunction reduction and P6 BLP postprocessing
%   post_bold - P7 BOLD postprocessing and P8 BLP-BOLD coupling
%   all       - all MATLAB-side phases above

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

local_disable_interactive_figures();

if ~exist('dry_run', 'var') || isempty(dry_run)
    dry_run = true;
end
if ~exist('backfill_phase', 'var') || isempty(backfill_phase)
    backfill_phase = 'all';
end
if ~exist('continue_on_error', 'var') || isempty(continue_on_error)
    continue_on_error = true;
end
if ~exist('minimal_bold_modes', 'var') || isempty(minimal_bold_modes)
    minimal_bold_modes = {'global_svd100', 'gsvd100_ds', 'HP_svd100', 'roi_mean'};
end
if ~exist('p3_bold_modes', 'var') || isempty(p3_bold_modes)
    p3_bold_modes = {'global_svd100', 'gsvd100_ds', 'HP_svd100', 'roi_mean'};
end
if ~exist('p5_methods', 'var') || isempty(p5_methods)
    p5_methods = {'svd', 'nmf', 'mds', 'umap'};
end
if ~exist('p5_component_count', 'var') || isempty(p5_component_count)
    p5_component_count = 3:8;
end
if ~exist('threshold_ratio', 'var') || isempty(threshold_ratio)
    threshold_ratio = 0.7;
end
if ~exist('backfill_datasets', 'var') || isempty(backfill_datasets)
    backfill_datasets = struct( ...
        'cfg_name', {'E10gb1', 'E10fV1', 'E10gH1', 'F12m01', 'E10gW1'}, ...
        'stem', {'e10gb1', 'e10fV1', 'e10gh1', 'f12m01', 'e10gw1'}, ...
        'dataset_id', {'E10.gb1', 'E10.fV1', 'E10.gH1', 'F12.m01', 'E10.gW1'});
end

dry_run = logical(dry_run);
continue_on_error = logical(continue_on_error);
backfill_phase = char(string(backfill_phase));
minimal_bold_modes = cellstr(string(minimal_bold_modes(:)).');
p3_bold_modes = cellstr(string(p3_bold_modes(:)).');
p5_methods = cellstr(string(p5_methods(:)).');
p5_component_count = double(p5_component_count);
threshold_ratio = double(threshold_ratio);

processed_root = io_project.get_project_processed_root();
autodl_blp_root = fullfile('E:', 'autodl_results_new');
autodl_bold_roots = {fullfile('E:', 'autodl_results_local', 'bold_wsl'), ...
    fullfile('E:', 'autodl_results', 'bold')};

fprintf('\nCurrent-five minimal backfill\n');
fprintf('  dry_run: %d\n', dry_run);
fprintf('  phase: %s\n', backfill_phase);
fprintf('  processed_root: %s\n', processed_root);
fprintf('  P3 BOLD modes: %s\n', strjoin(p3_bold_modes, ', '));
fprintf('  minimal BOLD modes: %s\n', strjoin(minimal_bold_modes, ', '));
fprintf('  P5 methods/k: %s / k=%s\n\n', ...
    strjoin(p5_methods, ', '), strjoin(cellstr(string(p5_component_count)), ', '));

phase_set = local_expand_phase(backfill_phase);
action_rows = table();

for i_dataset = 1:numel(backfill_datasets)
    item = backfill_datasets(i_dataset);
    fprintf('\n%s (%s)\n', item.dataset_id, item.stem);
    fprintf('%s\n', repmat('-', 1, 80));

    if any(strcmp(phase_set, 'pre_bold'))
        action_rows = [action_rows; local_backfill_pre_bold(item, processed_root, ...
            p3_bold_modes, dry_run, continue_on_error)]; %#ok<AGROW>
    end

    if any(strcmp(phase_set, 'blp'))
        action_rows = [action_rows; local_backfill_blp(item, processed_root, ...
            autodl_blp_root, p5_methods, p5_component_count, threshold_ratio, ...
            dry_run, continue_on_error)]; %#ok<AGROW>
    end

    if any(strcmp(phase_set, 'post_bold'))
        action_rows = [action_rows; local_backfill_post_bold(item, processed_root, ...
            autodl_bold_roots, minimal_bold_modes, dry_run, continue_on_error)]; %#ok<AGROW>
    end
end

fprintf('\nBackfill action summary\n');
fprintf('%s\n', repmat('=', 1, 80));
if isempty(action_rows)
    fprintf('No missing MATLAB-side actions were found for phase=%s.\n', backfill_phase);
else
    disp(action_rows);
end


function phases = local_expand_phase(phase)
switch lower(strtrim(phase))
    case 'all'
        phases = {'pre_bold', 'blp', 'post_bold'};
    case {'pre_bold', 'blp', 'post_bold'}
        phases = {lower(strtrim(phase))};
    otherwise
        error('Unknown backfill_phase: %s', phase);
end
end


function rows = local_backfill_pre_bold(item, processed_root, minimal_bold_modes, dry_run, continue_on_error)
rows = table();

missing_obs = local_missing_bold_observable_modes(processed_root, item, minimal_bold_modes);
if ~isempty(missing_obs)
    rows = [rows; local_action_row(item, 'P3', ...
        ['build BOLD observables: ' strjoin(missing_obs, ', ')])]; %#ok<AGROW>
    run_vars = struct();
    run_vars.cfg_name = item.cfg_name;
    run_vars.observable_modes = missing_obs;
    run_vars.force_recompute = false;
    local_run_entry('scripts/script_build_one_cfg_bold_observables.m', run_vars, dry_run, continue_on_error);
else
    fprintf('P3 minimal BOLD observables present; skipping build.\n');
end

missing_qc = local_missing_bold_qc_modes(processed_root, item, minimal_bold_modes);
if ~isempty(missing_qc)
    rows = [rows; local_action_row(item, 'P3_QC', ...
        ['plot BOLD QC: ' strjoin(missing_qc, ', ')])]; %#ok<AGROW>
    run_vars = struct();
    run_vars.cfg_name = item.cfg_name;
    run_vars.observable_modes = missing_qc;
    run_vars.save_qc = true;
    run_vars.qc_visible = 'off';
    local_run_entry('scripts/script_plot_one_cfg_bold_pre_reskoopnet_qc.m', run_vars, dry_run, continue_on_error);
else
    fprintf('P3 minimal BOLD QC present; skipping QC.\n');
end
end


function rows = local_backfill_blp(item, processed_root, autodl_blp_root, p5_methods, p5_component_count, threshold_ratio, dry_run, continue_on_error)
rows = table();

if ~local_has_p5_minimal(processed_root, item.stem, p5_methods, p5_component_count, threshold_ratio)
    rows = [rows; local_action_row(item, 'P5', 'run minimal consensus/density reduction')]; %#ok<AGROW>
    run_vars = struct();
    run_vars.cfg_name = item.cfg_name;
    run_vars.autodl_root = autodl_blp_root;
    run_vars.method_filter = p5_methods;
    run_vars.component_count_sweep = p5_component_count;
    run_vars.time_component_count_sweep = p5_component_count;
    run_vars.spectrum_component_count_sweep = p5_component_count;
    run_vars.make_overview_plot = false;
    run_vars.make_state_space_plot = true;
    run_vars.make_consensus_state_space_plot = true;
    run_vars.make_spectrum_diagnostics = true;
    run_vars.make_top30_window_plots = false;
    run_vars.make_thresholded_density = true;
    run_vars.make_thresholded_events = false;
    run_vars.make_dimred_thresholded_density = true;
    run_vars.make_dimred_thresholded_events = false;
    run_vars.threshold_mode = 'quantile';
    run_vars.threshold_ratio = threshold_ratio;
    run_vars.threshold_ratio_sweep = threshold_ratio;
    run_vars.continue_on_error = true;
    local_run_entry('scripts/script_run_one_cfg_blp_eigenfunction_reduction.m', run_vars, dry_run, continue_on_error);
else
    fprintf('P5 minimal consensus/density set present; skipping.\n');
end

if ~local_has_p6_outputs(processed_root, item.stem)
    rows = [rows; local_action_row(item, 'P6', 'run BLP eigenfunction postprocessing')]; %#ok<AGROW>
    run_vars = struct();
    run_vars.cfg_name = item.cfg_name;
    run_vars.autodl_root = autodl_blp_root;
    run_vars.force_recompute = false;
    run_vars.skip_existing = true;
    run_vars.max_basis = 30;
    run_vars.timescale_max_modes_sel = 30;
    run_vars.deconv_max_modes_sel = 30;
    run_vars.make_main_plot = true;
    run_vars.make_timescale_plot = true;
    run_vars.make_deconv_plot = true;
    run_vars.make_deconv_window_norm_plot = true;
    run_vars.make_spkt_residual_cross_correlation = true;
    run_vars.make_mua_residual_cross_correlation = true;
    local_run_entry('scripts/script_run_one_cfg_blp_eigenfunction_postprocessing.m', run_vars, dry_run, continue_on_error);
else
    fprintf('P6 outputs present; skipping.\n');
end
end


function rows = local_backfill_post_bold(item, processed_root, autodl_bold_roots, minimal_bold_modes, dry_run, continue_on_error)
rows = table();

available_modes = local_available_p4_bold_modes(autodl_bold_roots, item.stem, minimal_bold_modes);
if isempty(available_modes)
    fprintf('No P7-ready minimal P4 BOLD output chunks available; skipping P7/P8. Run or export P4 BOLD first.\n');
    return;
end

missing_p7 = local_missing_p7_modes(processed_root, item.stem, available_modes);
if ~isempty(missing_p7)
    rows = [rows; local_action_row(item, 'P7', ...
        ['run BOLD postprocessing: ' strjoin(missing_p7, ', ')])]; %#ok<AGROW>
    run_vars = struct();
    run_vars.cfg_name = item.cfg_name;
    run_vars.observable_modes = missing_p7;
    run_vars.force_recompute = false;
    run_vars.make_main_plot = true;
    run_vars.compute_deconv = true;
    run_vars.make_deconv_plot = true;
    run_vars.make_timescale_plot = true;
    run_vars.make_intrinsic_activation_maps = false;
    run_vars.make_intrinsic_roi_summary = false;
    run_vars.close_figures_after_each_run = true;
    local_run_entry('scripts/script_run_one_cfg_bold_reskoopnet_postprocessing.m', run_vars, dry_run, continue_on_error);
else
    fprintf('P7 minimal BOLD postprocessing present for available modes; skipping.\n');
end

missing_p8 = local_missing_p8_modes(processed_root, item.stem, available_modes);
if ~isempty(missing_p8)
    rows = [rows; local_action_row(item, 'P8', ...
        ['run BLP-BOLD coupling: ' strjoin(missing_p8, ', ')])]; %#ok<AGROW>
    run_vars = struct();
    run_vars.cfg_name = item.cfg_name;
    run_vars.observable_modes = missing_p8;
    run_vars.use_default_density_triplet = true;
    run_vars.require_all_density_sources = false;
    run_vars.output_mode = 'separate';
    run_vars.force_recompute = false;
    run_vars.max_lag_sec = 10;
    run_vars.border_mask_sec = 10;
    run_vars.top_n = 5;
    run_vars.make_xcorr_figures = false;
    run_vars.make_activation_maps = false;
    run_vars.make_roi_summaries = false;
    run_vars.force_activation_redraw = false;
    run_vars.activation_export_combined = true;
    run_vars.activation_export_by_density = true;
    run_vars.close_figures_after_each_run = true;
    local_run_entry('scripts/script_run_one_cfg_bold_cross_modal_coupling.m', run_vars, dry_run, continue_on_error);
else
    fprintf('P8 minimal BLP-BOLD coupling present for available modes; skipping.\n');
end
end


function missing = local_missing_bold_observable_modes(processed_root, item, modes)
missing = {};
stage_dir = fullfile(processed_root, item.stem, 'pipeline3_bold_observables');
for i = 1:numel(modes)
    file = fullfile(stage_dir, sprintf('%s_bold_observables_%s.mat', item.dataset_id, modes{i}));
    if exist(file, 'file') ~= 2
        missing{end + 1} = modes{i}; %#ok<AGROW>
    end
end
end


function missing = local_missing_bold_qc_modes(processed_root, item, modes)
missing = {};
stage_dir = fullfile(processed_root, item.stem, 'pipeline3_figures_bold_pre_reskoopnet_qc');
for i = 1:numel(modes)
    mode_dir = fullfile(stage_dir, modes{i});
    pngs = dir(fullfile(mode_dir, '*.png'));
    if isempty(pngs)
        missing{end + 1} = modes{i}; %#ok<AGROW>
    end
end
end


function tf = local_has_p5_minimal(processed_root, stem, methods, component_count, threshold_ratio)
condition_tags = {'abs_projected_vlambda', 'complex_split_projected_vlambda'};
component_count = double(component_count(:).');
for i_condition = 1:numel(condition_tags)
    condition_tag = condition_tags{i_condition};
    summary_dir = fullfile(processed_root, stem, 'pipeline5_summary_figures', ...
        condition_tag, 'eigenfunction_reduction', 'consensus_state_space');
    for i_count = 1:numel(component_count)
        k = component_count(i_count);
        for i = 1:numel(methods)
            pat = sprintf('%s_k%02d__efun_ssc__*.png', methods{i}, k);
            if isempty(dir(fullfile(summary_dir, pat)))
                tf = false;
                return;
            end
        end

        if ~local_has_p5_density(processed_root, stem, condition_tag, '', threshold_ratio)
            tf = false;
            return;
        end

        for i = 1:numel(methods)
            method_tag = sprintf('%s_k%02d', methods{i}, k);
            if ~local_has_p5_density(processed_root, stem, condition_tag, method_tag, threshold_ratio)
                tf = false;
                return;
            end
        end
    end
end
tf = true;
end


function tf = local_has_p5_density(processed_root, stem, condition_tag, method_tag, threshold_ratio)
q_tag = sprintf('q%03d', round(threshold_ratio * 100));
ratio_token = q_tag(2:end);
if isempty(method_tag)
    stage_dir = fullfile(processed_root, stem, 'pipeline5_raw_thresholded_density', ...
        condition_tag, 'mat');
    pat = sprintf('%s_*ratio_%s*%s*.mat', stem, ratio_token, condition_tag);
else
    stage_dir = fullfile(processed_root, stem, 'pipeline5_dimred_thresholded_density', ...
        condition_tag, method_tag, 'mat');
    pat = sprintf('%s_*ratio_%s*%s_%s*.mat', stem, ratio_token, condition_tag, method_tag);
end
tf = ~isempty(dir(fullfile(stage_dir, q_tag, pat))) || ...
    ~isempty(dir(fullfile(stage_dir, '**', pat)));
end


function tf = local_has_p6_outputs(processed_root, stem)
root = fullfile(processed_root, stem);
top_png = dir(fullfile(root, 'pipeline6_top_state_diversity_postprocessing', '**', '*.png'));
time_png = dir(fullfile(root, 'pipeline6_figures_timescale_diagnostics', '**', '*.png'));
spkt_mat = dir(fullfile(root, 'pipeline6_spkt_residual_cross_correlation', '**', '*.mat'));
mua_mat = dir(fullfile(root, 'pipeline6_mua_residual_cross_correlation', '**', '*.mat'));
tf = ~isempty(top_png) && ~isempty(time_png) && (~isempty(spkt_mat) || ~isempty(mua_mat));
end


function modes = local_available_p4_bold_modes(autodl_bold_roots, stem, requested_modes)
modes = {};
if ischar(autodl_bold_roots) || (isstring(autodl_bold_roots) && isscalar(autodl_bold_roots))
    autodl_bold_roots = {char(string(autodl_bold_roots))};
else
    autodl_bold_roots = cellstr(string(autodl_bold_roots(:)).');
end
for i = 1:numel(requested_modes)
    mode = requested_modes{i};
    found = false;
    for i_root = 1:numel(autodl_bold_roots)
        out_root = fullfile(autodl_bold_roots{i_root}, stem, 'mlp', 'outputs');
        if exist(out_root, 'dir') ~= 7
            continue;
        end
        run_dirs = dir(out_root);
        run_dirs = run_dirs([run_dirs.isdir]);
        for j = 1:numel(run_dirs)
            run_name = run_dirs(j).name;
            if startsWith(run_name, '.')
                continue;
            end
            if ~local_run_name_has_mode(run_name, mode)
                continue;
            end
            output_dir = fullfile(run_dirs(j).folder, run_name);
            if local_has_completed_p4_bold_run(output_dir)
                found = true;
                break;
            end
        end
        if found
            break;
        end
    end
    if found
        modes{end + 1} = mode; %#ok<AGROW>
    end
end
end


function tf = local_has_completed_p4_bold_run(output_dir)
summary_files = dir(fullfile(output_dir, '*summary.mat'));
if isempty(summary_files)
    tf = false;
    return;
end

chunk_files = dir(fullfile(output_dir, '*_outputs_*.mat'));
for i = 1:numel(chunk_files)
    name = chunk_files(i).name;
    if contains(name, '_outputs_Psi_')
        continue;
    end
    if ~isempty(regexp(name, '_outputs_\d+\.mat$', 'once'))
        tf = true;
        return;
    end
end
tf = false;
end


function missing = local_missing_p7_modes(processed_root, stem, modes)
missing = {};
p7_root = fullfile(processed_root, stem, 'pipeline7_bold_reskoopnet_postprocessing');
for i = 1:numel(modes)
    if ~local_has_p7_mode(p7_root, modes{i})
        missing{end + 1} = modes{i}; %#ok<AGROW>
    end
end
end


function tf = local_has_p7_mode(p7_root, mode)
tf = false;
if exist(p7_root, 'dir') ~= 7
    return;
end
dirs = dir(p7_root);
dirs = dirs([dirs.isdir]);
for i = 1:numel(dirs)
    name = dirs(i).name;
    if startsWith(name, '.') || ~local_run_name_has_mode(name, mode)
        continue;
    end
    mat_files = dir(fullfile(dirs(i).folder, name, '**', '*.mat'));
    png_files = dir(fullfile(dirs(i).folder, name, '**', '*.png'));
    if ~isempty(mat_files) && ~isempty(png_files)
        tf = true;
        return;
    end
end
end


function missing = local_missing_p8_modes(processed_root, stem, modes)
missing = {};
p8_root = fullfile(processed_root, stem, 'pipeline8_xcorr');
map_root = fullfile(processed_root, stem, 'pipeline8_top_maps');
for i = 1:numel(modes)
    tag = local_p8_tag_for_mode(modes{i});
    xcorr_dir = fullfile(p8_root, tag);
    map_dir = fullfile(map_root, tag);
    has_xcorr = ~isempty(dir(fullfile(xcorr_dir, '**', '*.csv'))) && ...
        ~isempty(dir(fullfile(xcorr_dir, '**', '*.mat')));
    has_maps = ~isempty(dir(fullfile(map_dir, '**', '*.png'))) && ...
        ~isempty(dir(fullfile(map_dir, '**', '*.mat')));
    if ~(has_xcorr && has_maps)
        missing{end + 1} = modes{i}; %#ok<AGROW>
    end
end
end


function tag = local_p8_tag_for_mode(mode)
switch lower(mode)
    case 'svd'
        tag = 'pv_svd';
    case 'global_svd100'
        tag = 'pv_gsvd100';
    case 'hp_svd100'
        tag = 'pv_hp100';
    case 'hp'
        tag = 'pv_hp';
    case 'elehp'
        tag = 'pv_elehp';
    case 'roi_mean'
        tag = 'pv_roi';
    case 'roi_mean_slow_band_power'
        tag = 'pv_roi_sbp';
    case 'slow_band_power_svd'
        tag = 'pv_sbp_svd';
    case 'gsvd100_ds'
        tag = 'pv_gsvd100_ds';
    otherwise
        tag = ['pv_' regexprep(lower(mode), '[^a-z0-9]+', '_')];
end
end


function tf = local_run_name_has_mode(run_name, mode)
suffix = ['projected_vlambda_' mode];
tf = endsWith(run_name, suffix);
end


function local_run_entry(script_rel, run_vars, dry_run, continue_on_error)
script_path = local_repo_path(script_rel);
msg = sprintf('MATLAB entry: %s', script_path);
local_disable_interactive_figures();
if dry_run
    fprintf('[dry-run] %s\n', msg);
    names = fieldnames(run_vars);
    for i = 1:numel(names)
        fprintf('          %s = %s\n', names{i}, local_value_preview(run_vars.(names{i})));
    end
    return;
end

try
    local_disable_interactive_figures();
    names = fieldnames(run_vars);
    for i = 1:numel(names)
        var_name = names{i}; %#ok<NASGU>
        var_value = run_vars.(names{i});
        eval([var_name ' = var_value;']);
    end
    run(script_path);
    local_disable_interactive_figures();
catch ME
    fprintf(2, 'Backfill entry failed: %s\n%s\n', script_path, getReport(ME, 'extended', 'hyperlinks', 'off'));
    if ~continue_on_error
        rethrow(ME);
    end
end
end


function local_disable_interactive_figures()
try
    set(groot, 'defaultFigureVisible', 'off');
    if isprop(groot, 'defaultAxesToolbarVisible')
        set(groot, 'defaultAxesToolbarVisible', 'off');
    end
catch
end
end


function path_out = local_repo_path(path_in)
if isfolder(fileparts(path_in)) && (startsWith(path_in, filesep) || ~isempty(regexp(path_in, '^[A-Za-z]:[\\/]', 'once')))
    path_out = path_in;
    return;
end
this_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(this_dir);
path_out = fullfile(repo_dir, path_in);
end


function row = local_action_row(item, stage, action)
row = table(string(item.dataset_id), string(item.stem), string(stage), string(action), ...
    'VariableNames', {'Dataset', 'Stem', 'Stage', 'Action'});
end


function text = local_value_preview(value)
if iscell(value)
    text = ['{' strjoin(cellstr(string(value(:)).'), ', ') '}'];
elseif ischar(value) || (isstring(value) && isscalar(value))
    text = char(string(value));
elseif isnumeric(value) || islogical(value)
    text = mat2str(value);
else
    text = class(value);
end
end
