this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));
set(groot, 'defaultFigureVisible', 'off');

cfg = cfg_eigenfunction_reduction_minimal();

fprintf('Eigenfunction reduction source:\n  %s\n', cfg.source.data_dir);
fprintf('Saving compact MAT result to:\n  %s\n', cfg.save.dir);
if isfield(cfg, 'viz') && isfield(cfg.viz, 'overview')
    fprintf('Figure save directory:\n  %s\n', cfg.viz.overview.save_dir);
end
fprintf('Save payload: %s\n', cfg.save.payload);

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

fig_state_space = [];
state_space_save_path = '';
if isfield(cfg, 'viz') && isfield(cfg.viz, 'state_space') && cfg.viz.state_space.enable
    [fig_state_space, state_space_info] = plot_eigenfunction_state_space_trajectory(result, cfg.viz.state_space);
    if isfield(state_space_info, 'save_path')
        state_space_save_path = state_space_info.save_path;
    end
end

fig_state_space_consensus = [];
state_space_consensus_save_path = '';
consensus_info = [];
source_consensus_file = '';
if isfield(cfg, 'viz') && isfield(cfg.viz, 'state_space_consensus') && ...
        cfg.viz.state_space_consensus.enable
    consensus_loader_cfg = struct();
    consensus_loader_cfg.file_stem = cfg.dataset.name;

    consensus_input = [];
    if isfield(cfg.viz.state_space_consensus, 'consensus_input')
        consensus_input = cfg.viz.state_space_consensus.consensus_input;
    end

    [C_consensus, source_consensus_file] = load_consensus_state_results( ...
        consensus_loader_cfg, cfg.dataset.processed_root, consensus_input);

    [fig_state_space_consensus, consensus_info] = ...
        plot_eigenfunction_state_space_consensus_trajectory( ...
            result, C_consensus, cfg.viz.state_space_consensus);

    if isfield(consensus_info, 'save_path')
        state_space_consensus_save_path = consensus_info.save_path;
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

print_artifact_paths(result.artifacts, 'thresholded_density_mat_file', ...
    'Saved thresholded density to:');
print_artifact_paths(result.artifacts, 'thresholded_density_figure_file', ...
    'Saved thresholded density figure to:');
print_artifact_paths(result.artifacts, 'thresholded_events_mat_file', ...
    'Saved thresholded events to:');
print_artifact_paths(result.artifacts, 'thresholded_events_figure_file', ...
    'Saved thresholded events figure to:');
print_artifact_paths(result.artifacts, 'dimred_thresholded_density_mat_file', ...
    'Saved dimension-reduced thresholded density to:');
print_artifact_paths(result.artifacts, 'dimred_thresholded_density_figure_file', ...
    'Saved dimension-reduced thresholded density figure to:');
print_artifact_paths(result.artifacts, 'dimred_thresholded_events_mat_file', ...
    'Saved dimension-reduced thresholded events to:');
print_artifact_paths(result.artifacts, 'dimred_thresholded_events_figure_file', ...
    'Saved dimension-reduced thresholded events figure to:');

if ~isempty(overview_save_path)
    fprintf('Saved overview figure to:\n  %s\n', overview_save_path);
end

if ~isempty(spectrum_save_path)
    fprintf('Saved spectrum diagnostics figure to:\n  %s\n', spectrum_save_path);
end

if ~isempty(state_space_save_path)
    fprintf('Saved state-space figure to:\n  %s\n', state_space_save_path);
end

if ~isempty(source_consensus_file)
    fprintf('Loaded consensus states from:\n  %s\n', source_consensus_file);
end

if ~isempty(state_space_consensus_save_path)
    fprintf('Saved consensus-state trajectory figure to:\n  %s\n', state_space_consensus_save_path);
end
if ~isempty(concat_info)
    fprintf('Concatenated %d chunks into %d samples.\n', ...
        concat_info.n_chunks, concat_info.total_length);
end

clear EDMD_outputs concat_info source_info fig_overview fig_spectrum fig_state_space fig_state_space_consensus
clear overview_info spectrum_info state_space_info consensus_info C_consensus source_consensus_file

function print_artifact_paths(artifacts, field_name, label)
if ~isfield(artifacts, field_name) || isempty(artifacts.(field_name))
    return;
end

paths = artifacts.(field_name);
if iscell(paths)
    paths = paths(:);
    keep = ~cellfun(@isempty, paths);
    paths = paths(keep);
    if isempty(paths)
        return;
    end
    fprintf('%s\n', label);
    for i = 1:numel(paths)
        fprintf('  %s\n', paths{i});
    end
else
    fprintf('%s\n  %s\n', label, paths);
end
end
