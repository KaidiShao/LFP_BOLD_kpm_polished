function summary_result = sync_one_pipeline8_cross_modal_summary_figures(run_info, source, output_root)
%SYNC_ONE_PIPELINE8_CROSS_MODAL_SUMMARY_FIGURES Copy key pipeline 8 PNGs into shared summaries.

if nargin < 1 || ~isstruct(run_info)
    error('run_info struct is required.');
end
if nargin < 2
    source = struct();
end
if nargin < 3 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

dataset_stem = char(string(local_get_field(run_info, 'dataset_stem', '')));
run_name = char(string(local_get_field(run_info, 'run_name', '')));
observable_mode = char(string(local_get_field(run_info, 'observable_mode', '')));
residual_form = char(string(local_get_field(run_info, 'residual_form', '')));
run_tag = char(string(local_get_field(run_info, 'run_tag', '')));
if isempty(run_tag)
    run_tag = make_bold_cross_modal_run_tag(run_info);
end
xcorr_dir = char(string(local_resolve_xcorr_dir(run_info, source)));

summary_result = struct();
summary_result.status = 'no_source_dir';
summary_result.scope = local_classify_observable_scope(observable_mode);
summary_result.dataset_stem = dataset_stem;
summary_result.run_name = run_name;
summary_result.run_tag = run_tag;
summary_result.observable_mode = observable_mode;
summary_result.residual_form = residual_form;
summary_result.summary_dir = '';
summary_result.copied_files = {};
summary_result.n_copied_files = 0;
summary_result.xcorr_dir = xcorr_dir;

if isempty(dataset_stem) || isempty(run_name) || isempty(observable_mode) || isempty(residual_form)
    summary_result.status = 'missing_run_info';
    return;
end
if isempty(xcorr_dir) || exist(xcorr_dir, 'dir') ~= 7
    return;
end

summary_root = io_project.get_pipeline_summary_dir(output_root, 8, 'figures_bold_xcorr_summary');
summary_dir = fullfile(summary_root, sprintf('%s_observables', summary_result.scope));
if exist(summary_dir, 'dir') ~= 7
    mkdir(summary_dir);
end

png_files = local_collect_pngs(xcorr_dir);
if isempty(png_files)
    summary_result.status = 'no_matching_pngs';
    summary_result.summary_dir = summary_dir;
    return;
end

copied_files = cell(numel(png_files), 1);
for i_file = 1:numel(png_files)
    src_file = png_files{i_file};
    [~, name, ext] = fileparts(src_file);
    dst_file = fullfile(summary_dir, sprintf('%s__%s%s', run_tag, name, ext));
    copyfile(src_file, dst_file, 'f');
    copied_files{i_file} = dst_file;
end

summary_result.status = 'ok';
summary_result.summary_dir = summary_dir;
summary_result.copied_files = copied_files;
summary_result.n_copied_files = numel(copied_files);
end


function xcorr_dir = local_resolve_xcorr_dir(run_info, source)
xcorr_dir = local_get_field(source, 'xcorr_dir', '');
if isempty(xcorr_dir)
    xcorr_dir = local_get_field(run_info, 'xcorr_dir', '');
end
if isempty(xcorr_dir)
    xcorr_file = local_get_field(source, 'xcorr_mat_file', '');
    if isempty(xcorr_file)
        xcorr_file = local_get_field(run_info, 'xcorr_file', '');
    end
    if ~isempty(xcorr_file)
        xcorr_dir = fileparts(char(string(xcorr_file)));
    end
end
xcorr_dir = char(string(xcorr_dir));
end


function scope = local_classify_observable_scope(observable_mode)
global_modes = {'svd', 'global_svd100', 'gsvd100_ds', 'global_slow_band_power_svd100', ...
    'roi_mean', 'roi_mean_slow_band_power'};
if any(strcmpi(observable_mode, global_modes))
    scope = 'global';
else
    scope = 'local';
end
end


function png_files = local_collect_pngs(xcorr_dir)
patterns = { ...
    fullfile(xcorr_dir, 'density', '*_summary.png'), ...
    fullfile(xcorr_dir, 'density', '*_top_signal_overlay.png'), ...
    fullfile(xcorr_dir, 'by_density', '*_summary.png'), ...
    fullfile(xcorr_dir, 'by_density', '*_top_signal_overlay.png')};
png_files = {};
for i_pat = 1:numel(patterns)
    L = dir(patterns{i_pat});
    for i_file = 1:numel(L)
        png_files{end + 1, 1} = fullfile(L(i_file).folder, L(i_file).name); %#ok<AGROW>
    end
end
if isempty(png_files)
    return;
end
png_files = unique(png_files, 'stable');
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
