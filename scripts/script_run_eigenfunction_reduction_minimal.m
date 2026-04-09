this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

cfg = cfg_eigenfunction_reduction_minimal();
[result, EDMD_outputs, concat_info, source_info] = run_eigenfunction_reduction_pipeline(cfg);

fig_overview = [];
overview_save_path = '';
if isfield(cfg, 'viz') && isfield(cfg.viz, 'overview') && cfg.viz.overview.enable
    [fig_overview, overview_info] = plot_eigenfunction_component_overview(result, cfg.viz.overview);
    if isfield(overview_info, 'save_path')
        overview_save_path = overview_info.save_path;
    end
end

fig_spectrum = [];
spectrum_save_path = '';
if strcmpi(result.meta.path_kind, 'spectrum') && ...
        isfield(cfg, 'viz') && isfield(cfg.viz, 'spectrum') && cfg.viz.spectrum.enable
    [fig_spectrum, spectrum_info] = plot_eigenfunction_spectrum_diagnostics(result, cfg.viz.spectrum);
    if isfield(spectrum_info, 'save_path')
        spectrum_save_path = spectrum_info.save_path;
    end
end

fprintf('Loaded EDMD source mode: %s\n', source_info.mode);
fprintf('Input size: T=%d, N=%d\n', ...
    size(result.data.efun_feature_time_by_mode, 1), ...
    size(result.data.efun_feature_time_by_mode, 2));
fprintf('Path kind: %s\n', result.meta.path_kind);
fprintf('Method: %s\n', result.meta.method);
fprintf('Temporal components size: [%s]\n', ...
    num2str(size(result.core.temporal_components_time_by_comp)));
fprintf('Mode weights size: [%s]\n', ...
    num2str(size(result.core.mode_weights_mode_by_comp)));

if isfield(result.artifacts, 'result_mat_file') && ~isempty(result.artifacts.result_mat_file)
    fprintf('Saved result to:\n  %s\n', result.artifacts.result_mat_file);
end

if ~isempty(overview_save_path)
    fprintf('Saved overview figure to:\n  %s\n', overview_save_path);
end

if ~isempty(spectrum_save_path)
    fprintf('Saved spectrum diagnostics figure to:\n  %s\n', spectrum_save_path);
end

if ~isempty(concat_info)
    fprintf('Concatenated %d chunks into %d samples.\n', ...
        concat_info.n_chunks, concat_info.total_length);
end

clear EDMD_outputs concat_info source_info fig_overview fig_spectrum overview_info spectrum_info
