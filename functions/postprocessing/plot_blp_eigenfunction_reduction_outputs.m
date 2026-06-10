function plot_paths = plot_blp_eigenfunction_reduction_outputs(result, cfg, params, C_consensus)
%PLOT_BLP_EIGENFUNCTION_REDUCTION_OUTPUTS Export canonical pipeline 5 figures for one method.

if nargin < 4
    C_consensus = [];
end

plot_paths = struct();
plot_paths.overview_png = "";
plot_paths.state_space_png = "";
plot_paths.consensus_state_space_png = "";
plot_paths.spectrum_diagnostics_png = "";

if cfg.viz.state_space.enable
    [fig, info] = plot_eigenfunction_state_space_trajectory(result, cfg.viz.state_space);
    if isfield(info, 'save_path')
        plot_paths.state_space_png = string(info.save_path);
    end
    local_close_if_requested(fig, params);
end

if cfg.viz.state_space_consensus.enable
    if isempty(C_consensus)
        error('Consensus-state trajectory plotting is enabled but no consensus input was provided.');
    end
    [fig, info] = plot_eigenfunction_state_space_consensus_trajectory( ...
        result, C_consensus, cfg.viz.state_space_consensus);
    if isfield(info, 'save_path')
        plot_paths.consensus_state_space_png = string(info.save_path);
    end
    local_close_if_requested(fig, params);
end

if strcmpi(result.meta.path_kind, 'spectrum') && cfg.viz.spectrum.enable
    [fig, info] = plot_eigenfunction_spectrum_diagnostics(result, cfg.viz.spectrum);
    if isfield(info, 'save_path')
        plot_paths.spectrum_diagnostics_png = string(info.save_path);
    end
    local_close_if_requested(fig, params);
end
end


function local_close_if_requested(fig, params)
if params.close_figures && ~isempty(fig) && isvalid(fig)
    close(fig);
end
end
