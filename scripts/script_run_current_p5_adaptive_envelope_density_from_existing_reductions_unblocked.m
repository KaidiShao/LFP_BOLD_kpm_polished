% Compute current P5 adaptive RMS-envelope density without re-running reduction.
%
% Raw density is recomputed from the current BLP EDMD outputs.  Dimred density
% reuses existing P5 reduction MAT files from the unsuffixed condition tags and
% writes new density MAT files under *_rmsenv_adaptive condition tags.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

if ~exist('dataset_filter', 'var'), dataset_filter = {}; end
if ~exist('condition_filter', 'var'), condition_filter = {}; end
if ~exist('method_filter', 'var') || isempty(method_filter)
    method_filter = {'svd', 'nmf', 'mds', 'umap'};
end
if ~exist('component_count_sweep', 'var') || isempty(component_count_sweep)
    component_count_sweep = 3:8;
end
if ~exist('activity_suffix', 'var') || isempty(activity_suffix)
    activity_suffix = 'rmsenv_adaptive';
end
if ~exist('source_condition_extra_suffix', 'var')
    source_condition_extra_suffix = '';
end
if ~exist('output_condition_extra_suffix', 'var') || isempty(output_condition_extra_suffix)
    output_condition_extra_suffix = activity_suffix;
end
if ~exist('autodl_root_override', 'var') || isempty(autodl_root_override)
    autodl_root_override = 'E:\autodl_results_new';
end
if ~exist('run_name_filter_override', 'var')
    run_name_filter_override = {};
end
if ~exist('dry_run', 'var') || isempty(dry_run)
    dry_run = false;
end
if ~exist('run_raw_density', 'var') || isempty(run_raw_density)
    run_raw_density = true;
end
if ~exist('run_dimred_density', 'var') || isempty(run_dimred_density)
    run_dimred_density = true;
end
if ~exist('skip_existing_density', 'var') || isempty(skip_existing_density)
    skip_existing_density = false;
end
if ~exist('allow_summary_only_outputs', 'var') || isempty(allow_summary_only_outputs)
    allow_summary_only_outputs = false;
end

if exist('dataset_specs_override', 'var') && ~isempty(dataset_specs_override)
    dataset_specs = dataset_specs_override;
else
    dataset_specs = struct( ...
        'dataset', {'e10gb1', 'e10fV1', 'e10gh1', 'e10gw1', 'f12m01', ...
            'k13m17', 'k13m23', 'f12m02', 'f12m03', 'f12m05'}, ...
        'cfg_name', {'E10gb1', 'E10fV1', 'E10gH1', 'E10gW1', 'F12m01', ...
            'K13m17', 'K13m23', 'F12m02', 'F12m03', 'F12m05'}, ...
        'conditions', {{'abs', 'complex_split'}, {'abs', 'complex_split'}, ...
            {'abs', 'complex_split'}, {'abs', 'complex_split'}, ...
            {'abs', 'complex_split'}, {'abs', 'complex_split'}, ...
            {'abs', 'complex_split'}, {'abs', 'complex_split'}, ...
            {'abs', 'complex_split'}, {'abs', 'complex_split'}} ...
        );
end

dataset_filter_lc = lower(string(cellstr(string(dataset_filter(:)).')));
condition_filter_lc = lower(string(cellstr(string(condition_filter(:)).')));

fprintf('P5 adaptive RMS-envelope density-only backfill\n');
fprintf('  reuse reductions: yes\n');
fprintf('  source suffix   : %s\n', char(string(source_condition_extra_suffix)));
fprintf('  output suffix   : %s\n', char(string(output_condition_extra_suffix)));
fprintf('  autodl root     : %s\n', char(string(autodl_root_override)));
fprintf('  methods         : %s\n', strjoin(cellstr(string(method_filter)), ', '));
fprintf('  k sweep         : %s\n', mat2str(component_count_sweep));
fprintf('  raw density     : %d\n', logical(run_raw_density));
fprintf('  dimred density  : %d\n', logical(run_dimred_density));
fprintf('  skip existing   : %d\n', logical(skip_existing_density));
fprintf('  summary-only    : %d\n', logical(allow_summary_only_outputs));
fprintf('  dry_run         : %d\n', logical(dry_run));

rows = struct([]);
for i_dataset = 1:numel(dataset_specs)
    spec = dataset_specs(i_dataset);
    dataset_lc = lower(string(spec.dataset));
    if ~isempty(dataset_filter_lc) && ~any(dataset_filter_lc == dataset_lc)
        continue;
    end

    cfg_name = spec.cfg_name;
    cfg = feval(['cfg_' cfg_name]);
    conditions = spec.conditions;
    if ~isempty(condition_filter_lc)
        keep = false(size(conditions));
        for i_cond = 1:numel(conditions)
            keep(i_cond) = any(condition_filter_lc == lower(string(conditions{i_cond})));
        end
        conditions = conditions(keep);
    end

    for i_cond = 1:numel(conditions)
        cond = conditions{i_cond};
        base_condition_tag = sprintf('%s_projected_vlambda', cond);
        source_condition_tag = local_append_suffix( ...
            base_condition_tag, source_condition_extra_suffix);
        output_condition_tag = local_append_suffix( ...
            base_condition_tag, output_condition_extra_suffix);
        fprintf('\n[P5 density-only] %s | %s -> %s\n', ...
            spec.dataset, source_condition_tag, output_condition_tag);

        params = local_build_density_params(spec.dataset, cond, output_condition_tag, ...
            method_filter, component_count_sweep, activity_suffix, ...
            autodl_root_override, run_name_filter_override, ...
            logical(allow_summary_only_outputs));
        runs = discover_completed_blp_mlp_runs(params);
        method_specs = build_blp_eigenfunction_reduction_method_specs(params);
        if isempty(runs)
            warning('No completed BLP MLP run found for %s | %s.', ...
                spec.dataset, cond);
            rows = local_append_row(rows, local_status_row( ...
                spec.dataset, output_condition_tag, '', 'discover_run', ...
                'missing_run', ''));
            continue;
        end
        if isempty(method_specs)
            warning('No P5 method specs found for %s | %s.', ...
                spec.dataset, cond);
            rows = local_append_row(rows, local_status_row( ...
                spec.dataset, output_condition_tag, '', 'method_specs', ...
                'missing_method_specs', ''));
            continue;
        end

        run_info = runs(1);
        if numel(runs) > 1
            fprintf('  [WARN] %d runs matched; using first: %s\n', ...
                numel(runs), run_info.run_name);
        end

        if ~run_raw_density
            rows = local_append_row(rows, local_status_row( ...
                spec.dataset, output_condition_tag, '', 'raw_density', ...
                'skipped', 'run_raw_density=false'));
        elseif skip_existing_density && local_has_raw_density( ...
                params.processed_root, spec.dataset, output_condition_tag)
            rows = local_append_row(rows, local_status_row( ...
                spec.dataset, output_condition_tag, '', 'raw_density', ...
                'skipped_existing', output_condition_tag));
        elseif dry_run
            rows = local_append_row(rows, local_status_row( ...
                spec.dataset, output_condition_tag, '', 'raw_density', ...
                'dry_run', run_info.output_dir));
        else
            try
                fprintf('  [raw] loading EDMD once for envelope density.\n');
                [EDMD_outputs, concat_info, source_info] = ...
                    load_blp_eigenfunction_reduction_source(run_info, params);
                raw_cfg = build_blp_raw_eigenfunction_threshold_cfg( ...
                    cfg, run_info, params, repo_root, ...
                    EDMD_outputs, concat_info, source_info);
                raw_cfg.thresholded_density.make_figure = false;
                raw_cfg.thresholded_density.save_figure = false;
                raw_cfg.thresholded_density.threshold_ratio_sweep = [];
                raw_cfg.thresholded_density.save_stem = sprintf('%s_td', ...
                    spec.dataset);
                raw_prep = prepare_blp_eigenfunction_reduction_inputs( ...
                    EDMD_outputs, raw_cfg);
                raw_result = build_blp_raw_eigenfunction_threshold_result( ...
                    raw_prep, run_info);
                run_blp_thresholded_density_stage( ...
                    raw_prep, raw_result, raw_cfg.thresholded_density);
                rows = local_append_row(rows, local_status_row( ...
                    spec.dataset, output_condition_tag, '', 'raw_density', ...
                    'ok', raw_cfg.thresholded_density.save_dir));
            catch ME_raw
                fprintf(2, '  [ERROR] raw density failed: %s\n', ME_raw.message);
                rows = local_append_row(rows, local_status_row( ...
                    spec.dataset, output_condition_tag, '', 'raw_density', ...
                    'error', ME_raw.message));
            end
            clear EDMD_outputs concat_info source_info raw_cfg raw_prep raw_result
        end

        if ~run_dimred_density
            rows = local_append_row(rows, local_status_row( ...
                spec.dataset, output_condition_tag, '', 'dimred_density', ...
                'skipped', 'run_dimred_density=false'));
            continue;
        end

        for i_method = 1:numel(method_specs)
            spec_method = method_specs(i_method);
            if ~spec_method.enabled
                continue;
            end
            [method_cfg, method_tag] = build_blp_eigenfunction_reduction_cfg( ...
                cfg, run_info, spec_method, params, repo_root, [], [], struct());
            if skip_existing_density && local_has_dimred_density( ...
                    params.processed_root, spec.dataset, ...
                    output_condition_tag, method_tag)
                rows = local_append_row(rows, local_status_row( ...
                    spec.dataset, output_condition_tag, method_tag, ...
                    'dimred_density', 'skipped_existing', ''));
                continue;
            end
            source_result_file = local_latest_reduction_result( ...
                params.processed_root, spec.dataset, source_condition_tag, method_tag);
            if isempty(source_result_file)
                fprintf(2, '  [WARN] missing source reduction: %s | %s\n', ...
                    source_condition_tag, method_tag);
                rows = local_append_row(rows, local_status_row( ...
                    spec.dataset, output_condition_tag, method_tag, ...
                    'dimred_density', 'missing_source_reduction', ''));
                continue;
            end

            if dry_run
                rows = local_append_row(rows, local_status_row( ...
                    spec.dataset, output_condition_tag, method_tag, ...
                    'dimred_density', 'dry_run', source_result_file));
                continue;
            end

            try
                fprintf('  [dimred] %s <- %s\n', method_tag, source_result_file);
                S = load(source_result_file, 'result');
                result = S.result;
                td_cfg = method_cfg.dimred_thresholded_density;
                td_cfg.make_figure = false;
                td_cfg.save_figure = false;
                td_cfg.save_results = true;
                td_cfg.threshold_ratio_sweep = [];
                % Keep the dimred density filename short.  Long standardized
                % condition tags plus umap_kXX can exceed Windows/MATLAB HDF5
                % path limits when using -v7.3.
                td_cfg.save_stem = sprintf('%s', spec.dataset);
                dt = [];
                if isfield(result, 'input') && isstruct(result.input) && ...
                        isfield(result.input, 'dt')
                    dt = result.input.dt;
                end
                get_dimred_thresholded_density(result, dt, td_cfg);
                rows = local_append_row(rows, local_status_row( ...
                    spec.dataset, output_condition_tag, method_tag, ...
                    'dimred_density', 'ok', td_cfg.save_dir));
            catch ME_dim
                fprintf(2, '  [ERROR] dimred density failed for %s: %s\n', ...
                    method_tag, ME_dim.message);
                rows = local_append_row(rows, local_status_row( ...
                    spec.dataset, output_condition_tag, method_tag, ...
                    'dimred_density', 'error', ME_dim.message));
            end
            clear method_cfg result S td_cfg
        end
    end
end

manifest_dir = fullfile(repo_root, 'results', 'manifests');
if exist(manifest_dir, 'dir') ~= 7
    mkdir(manifest_dir);
end
manifest_file = fullfile(manifest_dir, ...
    ['p5_rmsenv_adaptive_density_only_', char(datetime('now', ...
    'Format', 'yyyyMMdd_HHmmss')), '.csv']);
if ~isempty(rows)
    writetable(struct2table(rows(:), 'AsArray', true), manifest_file);
end
fprintf('\nFinished P5 adaptive RMS-envelope density-only backfill.\n');
fprintf('Manifest: %s\n', manifest_file);


function params = local_build_density_params(dataset, cond, output_condition_tag, ...
        method_filter, component_count_sweep, activity_suffix, ...
        autodl_root_override, run_name_filter_override, allow_summary_only_outputs)
params = build_blp_eigenfunction_reduction_params();
params.dataset_stems = {dataset};
params.autodl_root = char(string(autodl_root_override));
params.allow_summary_only_outputs = logical(allow_summary_only_outputs);
params.condition_key_filter = {sprintf('%s|%s|projected_vlambda', dataset, cond)};
if ~isempty(run_name_filter_override)
    params.run_name_filter = cellstr(string(run_name_filter_override(:)).');
end
params.condition_output_tag = output_condition_tag;
params.condition_tag_mode = 'condition';
params.method_filter = method_filter;
params.component_count_sweep = component_count_sweep;
params.do_thresholded_density = true;
params.do_thresholded_events = false;
params.do_dimred_thresholded_density = true;
params.do_dimred_thresholded_events = false;
params.threshold_mode = 'quantile';
params.threshold_ratio = 0.7;
params.threshold_ratio_sweep = [];
params.density_value_transform = 'abs';
params.lfp_activity_transform = 'rms_envelope';
params.lfp_activity_window_policy = 'eigenvalue_adaptive_rms_envelope';
params.envelope_enable = true;
params.envelope_policy = 'eigenvalue_adaptive_rms';
params.envelope_alpha = 0.35;
params.envelope_min_window_sec = 0.03;
params.envelope_max_window_sec = 1.0;
params.envelope_fallback_window_sec = 0.10;
params.component_timescale_weight_transform = 'abs';
params.component_timescale_top_modes = 5;
params.make_overview_plot = false;
params.make_state_space_plot = false;
params.make_consensus_state_space_plot = false;
params.make_spectrum_diagnostics = false;
params.make_top30_window_plots = false;
params.save_thresholded_density_results = true;
params.save_dimred_thresholded_density_results = true;
params.continue_on_error = true;
params.activity_suffix = activity_suffix;
end


function tag = local_append_suffix(base_tag, suffix)
tag = char(string(base_tag));
suffix = char(string(suffix));
if isempty(suffix)
    return;
end
tag = sprintf('%s_%s', tag, suffix);
tag = regexprep(tag, '[^a-zA-Z0-9_\-]+', '_');
tag = regexprep(tag, '_+', '_');
tag = regexprep(tag, '^_+|_+$', '');
end


function tf = local_has_raw_density(processed_root, dataset, condition_tag)
stage_dir = fullfile(processed_root, dataset, ...
    'pipeline5_raw_thresholded_density', condition_tag, 'mat');
pat = sprintf('%s_*ratio_070*%s*.mat', dataset, condition_tag);
tf = ~isempty(dir(fullfile(stage_dir, '**', pat)));
end


function tf = local_has_dimred_density(processed_root, dataset, condition_tag, method_tag)
stage_dir = fullfile(processed_root, dataset, ...
    'pipeline5_dimred_thresholded_density', condition_tag, method_tag, 'mat');
pat = sprintf('%s_*ratio_070*%s_%s*.mat', dataset, condition_tag, method_tag);
tf = ~isempty(dir(fullfile(stage_dir, '**', pat)));
end


function file_name = local_latest_reduction_result(processed_root, dataset, condition_tag, method_tag)
mat_dir = fullfile(processed_root, dataset, 'pipeline5_eigenfunction_reduction', ...
    condition_tag, method_tag, 'mat');
L = dir(fullfile(mat_dir, '*.mat'));
if isempty(L)
    file_name = '';
    return;
end
[~, idx] = max([L.datenum]);
file_name = fullfile(L(idx).folder, L(idx).name);
end


function rows = local_append_row(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end + 1) = row;
end
end


function row = local_status_row(dataset, condition_tag, method_tag, stage, status, detail)
row = struct();
row.dataset = string(dataset);
row.condition_tag = string(condition_tag);
row.method_tag = string(method_tag);
row.stage = string(stage);
row.status = string(status);
row.detail = string(detail);
row.created_at = string(datetime('now'));
end
