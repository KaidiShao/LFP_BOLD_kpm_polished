function hfig = plot_bold_segment(D, time_range_sec, params)
% Plot a BOLD segment from a loaded or preprocessed BOLD struct.

if nargin < 2 || isempty(time_range_sec)
    error('time_range_sec must be [t_start, t_end].');
end
if nargin < 3
    params = struct();
end

if ~isfield(params, 'source') || isempty(params.source)
    params.source = 'data';
end
if ~isfield(params, 'max_variables') || isempty(params.max_variables)
    params.max_variables = 20;
end
if ~isfield(params, 'trace_scale') || isempty(params.trace_scale)
    params.trace_scale = 1;
end

X = select_plot_source(D, params.source);
[t, ~, ~, ~] = build_global_time_axis_from_sessions(D.session_lengths, D.session_dx);

idx = find(t >= time_range_sec(1) & t <= time_range_sec(2));
if isempty(idx)
    error('Requested time range does not overlap BOLD data.');
end

n_plot = min(params.max_variables, size(X, 2));
x_seg = double(X(idx, 1:n_plot));
t_seg = t(idx);

x_disp = zeros(size(x_seg));
for j = 1:n_plot
    xj = x_seg(:, j);
    xj = xj - mean(xj, 'omitnan');
    sj = std(xj, 0, 'omitnan');
    if ~isfinite(sj) || sj == 0
        sj = 1;
    end
    x_disp(:, j) = params.trace_scale * xj / sj;
end

offsets = (n_plot-1:-1:0);

hfig = figure('Color', 'w', 'Position', [100, 100, 1100, 460]);
ax = axes('Parent', hfig);
hold(ax, 'on');

for j = 1:n_plot
    plot(ax, t_seg, x_disp(:, j) + offsets(j), 'k', 'LineWidth', 0.8);
end

border_t = t(D.border_idx);
border_t = border_t(border_t >= time_range_sec(1) & border_t <= time_range_sec(2));
for i = 1:numel(border_t)
    xline(ax, border_t(i), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end

if isfield(D, 'variable_labels') && numel(D.variable_labels) >= n_plot
    ylabels = D.variable_labels(1:n_plot);
else
    ylabels = cellstr(compose("var%04d", (1:n_plot).'));
end

set(ax, 'YTick', offsets, 'YTickLabel', ylabels);
set(ax, 'Box', 'off');
xlabel(ax, 'Time (s)');
title(ax, sprintf('BOLD segment: %.3f-%.3f s', time_range_sec(1), time_range_sec(2)));
xlim(ax, time_range_sec);
drawnow;
end

function X = select_plot_source(D, source)
switch lower(source)
    case 'data'
        X = D.data;
    case 'filtered'
        if isfield(D, 'filtered_data') && ~isempty(D.filtered_data)
            X = D.filtered_data;
        else
            X = D.data;
        end
    otherwise
        error('Unsupported plot source: %s', source);
end
end
