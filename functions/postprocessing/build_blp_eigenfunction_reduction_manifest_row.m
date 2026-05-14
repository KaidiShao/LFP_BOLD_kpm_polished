function row = build_blp_eigenfunction_reduction_manifest_row( ...
        run_info, method_tag, status, message, result, plot_paths, runtime_sec)
%BUILD_BLP_EIGENFUNCTION_REDUCTION_MANIFEST_ROW Summarize one reduction run.

row = struct();
row.dataset = string(run_info.dataset);
row.observable_mode = string(run_info.observable_mode);
row.residual_form = string(run_info.residual_form);
row.run_name = string(run_info.run_name);
row.method = string(method_tag);
row.status = string(status);
row.message = string(message);
row.runtime_sec = double(runtime_sec);
row.n_components = NaN;
row.result_mat_file = "";
row.overview_png = "";
row.state_space_png = "";
row.consensus_state_space_png = "";
row.spectrum_diagnostics_png = "";
row.top30_window_manifest_file = "";
row.reconstruction_error_fro = NaN;

if nargin < 5 || isempty(result)
    return;
end

if isfield(result, 'artifacts') && isfield(result.artifacts, 'result_mat_file')
    row.result_mat_file = string(result.artifacts.result_mat_file);
end
if isfield(result, 'core') && isfield(result.core, 'n_components')
    row.n_components = double(result.core.n_components);
end
if nargin >= 6 && ~isempty(plot_paths)
    if isfield(plot_paths, 'overview_png')
        row.overview_png = string(plot_paths.overview_png);
    end
    if isfield(plot_paths, 'state_space_png')
        row.state_space_png = string(plot_paths.state_space_png);
    end
    if isfield(plot_paths, 'consensus_state_space_png')
        row.consensus_state_space_png = string(plot_paths.consensus_state_space_png);
    end
    if isfield(plot_paths, 'spectrum_diagnostics_png')
        row.spectrum_diagnostics_png = string(plot_paths.spectrum_diagnostics_png);
    end
    if isfield(plot_paths, 'top30_window_manifest_file')
        row.top30_window_manifest_file = string(plot_paths.top30_window_manifest_file);
    end
end
if isfield(result, 'quality') && isfield(result.quality, 'reconstruction_error_fro')
    row.reconstruction_error_fro = double(result.quality.reconstruction_error_fro);
end
end
