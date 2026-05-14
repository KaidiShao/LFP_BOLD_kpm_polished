function [te_summary, artifacts] = run_blp_thresholded_events_stage(prep, result, te_cfg)
%RUN_BLP_THRESHOLDED_EVENTS_STAGE Run raw eigenfunction thresholded event stage.

if isempty(prep.dt)
    error('Thresholded events require a sampling interval. %s', ...
        prep.dt_source.message);
end

te_cfg.time_axis = prep.time_axis;
te_cfg.dt_source = prep.dt_source;
te_cfg.mode_index = result.input.mode_index;
te_cfg.selected_mode_idx_in_original = prep.selected_mode_idx_in_original;

if ~isfield(te_cfg, 'title') || isempty(te_cfg.title)
    te_cfg.title = sprintf('%s Thresholded Eigenfunction Events', ...
        result.meta.method);
end

[E, figs] = get_thresholded_events( ...
    prep.efun_feature_time_by_mode, prep.dt, prep.evalues_discrete, te_cfg);

if te_cfg.close_figure && isfield(figs, 'summary') && ...
        ~isempty(figs.summary) && isvalid(figs.summary)
    close(figs.summary);
end

te_summary = local_compact_thresholded_events(E);
artifacts = E.artifacts;
end


function te_summary = local_compact_thresholded_events(E)
te_summary = struct();
te_summary.created_at = E.created_at;
te_summary.input = E.input;
te_summary.params = E.params;
te_summary.summary = E.summary;
te_summary.artifacts = E.artifacts;
te_summary.threshold_by_mode = E.threshold_by_mode;
te_summary.mode_timescales = E.mode_timescales;
te_summary.event_rate_size = size(E.event_rate_time_by_mode);
te_summary.t_range = [E.t_centers(1), E.t_centers(end)];
te_summary.n_windows = numel(E.t_centers);
te_summary.n_modes = size(E.event_rate_time_by_mode, 2);
te_summary.n_events = height(E.event_table);
end
