function [te_summary, artifacts] = run_blp_dimred_thresholded_events_stage( ...
        prep, result, te_cfg)
%RUN_BLP_DIMRED_THRESHOLDED_EVENTS_STAGE Run thresholded event extraction on reduced components.

if isempty(result.core.temporal_components_time_by_comp)
    error('Dimension-reduced thresholded events require temporal components.');
end

if isempty(prep.dt)
    error(['Dimension-reduced thresholded events require a sampling ', ...
        'interval. %s'], prep.dt_source.message);
end

te_cfg.time_axis = prep.time_axis;
te_cfg.dt_source = prep.dt_source;
te_cfg.component_index = (1:size(result.core.temporal_components_time_by_comp, 2)).';

if ~isfield(te_cfg, 'title') || isempty(te_cfg.title)
    te_cfg.title = sprintf('%s Thresholded Component Events', ...
        result.meta.method);
end

[E, ~] = get_dimred_thresholded_events(result, prep.dt, te_cfg);

te_summary = local_compact_dimred_thresholded_events(E);
artifacts = E.artifacts;
end


function te_summary = local_compact_dimred_thresholded_events(E)
te_summary = struct();
te_summary.created_at = E.created_at;
te_summary.meta = E.meta;
te_summary.input = E.input;
te_summary.params = E.params;
te_summary.summary = E.summary;
te_summary.artifacts = E.artifacts;
te_summary.threshold_by_component = E.threshold_by_component;
te_summary.component_index = E.component_index;
te_summary.component_evalues_discrete = E.component_evalues_discrete;
te_summary.component_evalue_info = E.component_evalue_info;
te_summary.component_timescales = E.component_timescales;
te_summary.event_rate_size = size(E.event_rate_time_by_component);
te_summary.t_range = [E.t_centers(1), E.t_centers(end)];
te_summary.n_windows = numel(E.t_centers);
te_summary.n_components = size(E.event_rate_time_by_component, 2);
te_summary.n_events = height(E.event_table);
end
