function [td_summary, artifacts] = run_blp_thresholded_density_stage(prep, result, td_cfg)
%RUN_BLP_THRESHOLDED_DENSITY_STAGE Run raw eigenfunction thresholded density output stage.

if isempty(prep.dt) && local_thresholded_density_needs_dt(td_cfg)
    error(['Thresholded density requires a sampling interval because its ', ...
        'window/step/smoothing parameters are specified in seconds. ', ...
        '%s'], prep.dt_source.message);
end

td_cfg.time_axis = prep.time_axis;
td_cfg.dt_source = prep.dt_source;
td_cfg.mode_index = result.input.mode_index;
td_cfg.selected_mode_idx_in_original = prep.selected_mode_idx_in_original;

if ~isfield(td_cfg, 'title') || isempty(td_cfg.title)
    td_cfg.title = sprintf('%s Thresholded Eigenfunction Density', ...
        result.meta.method);
end

[td_cfgs, ratio_tags] = expand_blp_threshold_sweep_configs(td_cfg);
td_summary = struct([]);
artifact_rows = struct([]);

for i_cfg = 1:numel(td_cfgs)
    cfg_i = td_cfgs(i_cfg);
    [D, fig] = get_thresholded_density( ...
        prep.efun_feature_time_by_mode, prep.dt, cfg_i);

    if cfg_i.close_figure && ~isempty(fig) && isvalid(fig)
        close(fig);
    end

    row_i = local_compact_thresholded_density(D);
    row_i.threshold_sweep_tag = char(ratio_tags(i_cfg));
    row_i.threshold_ratio = double(cfg_i.threshold_ratio);
    td_summary = local_append_struct(td_summary, row_i);

    artifacts_i = D.artifacts;
    artifacts_i.threshold_sweep_tag = char(ratio_tags(i_cfg));
    artifacts_i.threshold_ratio = double(cfg_i.threshold_ratio);
    artifact_rows = local_append_struct(artifact_rows, artifacts_i);
end

artifacts = local_collect_sweep_artifacts(artifact_rows, td_cfgs, td_cfg);
end


function tf = local_thresholded_density_needs_dt(td_cfg)
has_window_samples = isfield(td_cfg, 'window_samples') && ~isempty(td_cfg.window_samples);
has_step_samples = isfield(td_cfg, 'step_samples') && ~isempty(td_cfg.step_samples);

tf = ~(has_window_samples && has_step_samples);

smooth_enabled = isfield(td_cfg, 'smooth_density') && ~isempty(td_cfg.smooth_density) && ...
    logical(td_cfg.smooth_density);
has_smooth_bins = isfield(td_cfg, 'smooth_window_bins') && ~isempty(td_cfg.smooth_window_bins);
smooth_sec_positive = isfield(td_cfg, 'smooth_window_sec') && ...
    ~isempty(td_cfg.smooth_window_sec) && td_cfg.smooth_window_sec > 0;

tf = tf || (smooth_enabled && smooth_sec_positive && ~has_smooth_bins);
end


function td_summary = local_compact_thresholded_density(D)
td_summary = struct();
td_summary.created_at = D.created_at;
td_summary.input = D.input;
td_summary.params = D.params;
td_summary.summary = D.summary;
td_summary.artifacts = D.artifacts;
td_summary.threshold_by_mode = D.threshold_by_mode;
td_summary.density_size = size(D.density_time_by_mode);
td_summary.t_range = [D.t_centers(1), D.t_centers(end)];
td_summary.n_windows = numel(D.t_centers);
td_summary.n_modes = size(D.density_time_by_mode, 2);
end


function rows = local_append_struct(rows, row_i)
if isempty(rows)
    rows = row_i;
else
    rows(end + 1) = row_i;
end
end


function artifacts = local_collect_sweep_artifacts(artifact_rows, cfgs, original_cfg)
if isempty(artifact_rows)
    artifacts = struct('mat_file', '', 'figure_file', '', ...
        'mat_files', strings(0, 1), 'figure_files', strings(0, 1), ...
        'threshold_sweep_tags', strings(0, 1), 'threshold_ratios', []);
    return;
end

primary_idx = local_primary_sweep_index(cfgs, original_cfg);
artifacts = artifact_rows(primary_idx);

n = numel(artifact_rows);
mat_files = strings(n, 1);
figure_files = strings(n, 1);
tags = strings(n, 1);
ratios = nan(n, 1);
for i = 1:n
    mat_files(i) = string(local_get_artifact_field(artifact_rows(i), 'mat_file'));
    figure_files(i) = string(local_get_artifact_field(artifact_rows(i), 'figure_file'));
    tags(i) = string(local_get_artifact_field(artifact_rows(i), 'threshold_sweep_tag'));
    ratios(i) = double(local_get_artifact_field(artifact_rows(i), 'threshold_ratio'));
end

artifacts.mat_files = mat_files;
artifacts.figure_files = figure_files;
artifacts.threshold_sweep_tags = tags;
artifacts.threshold_ratios = ratios;
end


function idx = local_primary_sweep_index(cfgs, original_cfg)
idx = 1;
if ~isfield(original_cfg, 'threshold_ratio') || isempty(original_cfg.threshold_ratio)
    return;
end

target = double(original_cfg.threshold_ratio);
for i = 1:numel(cfgs)
    if isfield(cfgs(i), 'threshold_ratio') && ...
            abs(double(cfgs(i).threshold_ratio) - target) < 1e-12
        idx = i;
        return;
    end
end
end


function value = local_get_artifact_field(S, field_name)
if isfield(S, field_name) && ~isempty(S.(field_name))
    value = S.(field_name);
else
    value = '';
end
end
