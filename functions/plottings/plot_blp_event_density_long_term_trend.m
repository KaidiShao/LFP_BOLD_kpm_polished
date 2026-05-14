function [hfig, info] = plot_blp_event_density_long_term_trend(input_data, params)
%PLOT_BLP_EVENT_DENSITY_LONG_TERM_TREND Plot all pipeline2 band densities in one figure.
%
%   [hfig, info] = plot_blp_event_density_long_term_trend(input_data, params)
%
%   input_data can be:
%     1) a saved pipeline2 event-density MAT file containing variable E
%     2) a struct E returned by compute_blp_event_density
%
%   The plot is intended as a quick long-timescale diagnostic to compare the
%   three canonical band-density traces within one dataset.

if nargin < 2
    params = struct();
end
params = local_fill_default_params(params);

[E, info] = local_load_event_density_input(input_data);
[density_to_plot, raw_density] = local_resolve_density_series(E, params);
t_centers = local_resolve_time_axis(E, size(density_to_plot, 1));
density_to_plot = local_apply_trend_smoothing(density_to_plot, E, t_centers, params);
[t_plot, time_unit_label, time_scale] = local_convert_time_axis(t_centers, params.time_unit);
[boundary_t, session_mid_t, session_labels] = local_resolve_session_time_info(E, time_scale);
band_labels = local_resolve_band_labels(E, size(density_to_plot, 2));
dataset_label = local_resolve_dataset_label(E, info);
info.dataset_id = local_get_optional_string_field(E, 'dataset_id');
info.file_stem = local_get_optional_string_field(E, 'file_stem');
info.total_duration_sec = local_get_total_duration_sec(E, t_centers);
info.trend_window_sec = double(params.trend_window_sec);

hfig = figure( ...
    'Color', 'w', ...
    'Visible', params.visible, ...
    'Position', params.figure_position);
ax = axes('Parent', hfig);
if exist('disableDefaultInteractivity', 'file') == 2
    disableDefaultInteractivity(ax);
end
hold(ax, 'on');

if params.show_session_boundaries
    for i = 1:numel(boundary_t)
        xline(ax, boundary_t(i), params.boundary_line_style, ...
            'Color', params.boundary_color, ...
            'LineWidth', params.boundary_line_width, ...
            'HandleVisibility', 'off');
    end
end

if params.show_raw_mean && ~isempty(raw_density)
    for b = 1:size(raw_density, 2)
        plot(ax, t_plot, raw_density(:, b), ...
            'Color', local_lighten_color(params.band_colors(b, :), params.raw_color_mix), ...
            'LineWidth', params.raw_line_width, ...
            'HandleVisibility', 'off');
    end
end

plot_handles = gobjects(size(density_to_plot, 2), 1);
for b = 1:size(density_to_plot, 2)
    plot_handles(b) = plot(ax, t_plot, density_to_plot(:, b), ...
        'Color', params.band_colors(b, :), ...
        'LineWidth', params.line_width);
end

if params.show_session_labels && ~isempty(session_mid_t)
    y_limits = ylim(ax);
    y_text = y_limits(2) - 0.03 * max(eps, diff(y_limits));
    for i = 1:numel(session_mid_t)
        text(ax, session_mid_t(i), y_text, session_labels{i}, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontSize', params.session_label_font_size, ...
            'Color', params.session_label_color, ...
            'Clipping', 'on');
    end
end

legend(ax, plot_handles, band_labels, 'Location', params.legend_location);
title(ax, local_build_title(dataset_label, E, info.total_duration_sec, params), 'Interpreter', 'none');
xlabel(ax, sprintf('Time (%s)', time_unit_label));
ylabel(ax, params.y_label);
xlim(ax, [t_plot(1), t_plot(end)]);
grid(ax, 'on');
set(ax, 'Box', 'off');

info.band_labels = string(band_labels(:));
info.boundary_t = boundary_t(:);
info.time_unit = string(time_unit_label);
info.time_scale = time_scale;
info.n_sessions = numel(session_labels);
info.output_files = strings(0, 1);

if params.save
    if exist(params.output_dir, 'dir') ~= 7
        mkdir(params.output_dir);
    end

    output_base = fullfile(params.output_dir, local_build_figure_prefix(E, info, params));
    info.output_files = local_save_figure(hfig, output_base, params.save_formats, params.resolution);
end

if params.close_after_save && params.save
    close(hfig);
end
end


function params = local_fill_default_params(params)
if ~isfield(params, 'save') || isempty(params.save)
    params.save = false;
end
if ~isfield(params, 'output_dir') || isempty(params.output_dir)
    params.output_dir = fullfile(pwd, 'blp_event_density_long_term_trends');
end
if ~isfield(params, 'figure_prefix')
    params.figure_prefix = '';
end
if ~isfield(params, 'visible') || isempty(params.visible)
    params.visible = 'on';
end
if ~isfield(params, 'close_after_save') || isempty(params.close_after_save)
    params.close_after_save = false;
end
if ~isfield(params, 'save_formats') || isempty(params.save_formats)
    params.save_formats = {'png'};
end
if ~isfield(params, 'resolution') || isempty(params.resolution)
    params.resolution = 220;
end
if ~isfield(params, 'use_smoothed') || isempty(params.use_smoothed)
    params.use_smoothed = true;
end
if ~isfield(params, 'show_raw_mean') || isempty(params.show_raw_mean)
    params.show_raw_mean = false;
end
if ~isfield(params, 'trend_window_sec') || isempty(params.trend_window_sec)
    params.trend_window_sec = 0;
end
if ~isfield(params, 'time_unit') || isempty(params.time_unit)
    params.time_unit = 'auto';
end
if ~isfield(params, 'show_session_boundaries') || isempty(params.show_session_boundaries)
    params.show_session_boundaries = true;
end
if ~isfield(params, 'show_session_labels') || isempty(params.show_session_labels)
    params.show_session_labels = false;
end
if ~isfield(params, 'band_colors') || isempty(params.band_colors)
    params.band_colors = [ ...
        0.0000, 0.4470, 0.7410; ...
        0.8500, 0.3250, 0.0980; ...
        0.4660, 0.6740, 0.1880];
end
if ~isfield(params, 'line_width') || isempty(params.line_width)
    params.line_width = 1.8;
end
if ~isfield(params, 'raw_line_width') || isempty(params.raw_line_width)
    params.raw_line_width = 0.8;
end
if ~isfield(params, 'raw_color_mix') || isempty(params.raw_color_mix)
    params.raw_color_mix = 0.65;
end
if ~isfield(params, 'boundary_color') || isempty(params.boundary_color)
    params.boundary_color = [0.78, 0.78, 0.78];
end
if ~isfield(params, 'boundary_line_style') || isempty(params.boundary_line_style)
    params.boundary_line_style = '--';
end
if ~isfield(params, 'boundary_line_width') || isempty(params.boundary_line_width)
    params.boundary_line_width = 0.9;
end
if ~isfield(params, 'legend_location') || isempty(params.legend_location)
    params.legend_location = 'eastoutside';
end
if ~isfield(params, 'y_label') || isempty(params.y_label)
    params.y_label = 'Event density (events/s)';
end
if ~isfield(params, 'title_suffix') || isempty(params.title_suffix)
    params.title_suffix = '';
end
if ~isfield(params, 'figure_position') || isempty(params.figure_position)
    params.figure_position = [100, 100, 1300, 460];
end
if ~isfield(params, 'session_label_font_size') || isempty(params.session_label_font_size)
    params.session_label_font_size = 8;
end
if ~isfield(params, 'session_label_color') || isempty(params.session_label_color)
    params.session_label_color = [0.35, 0.35, 0.35];
end
if ~isscalar(params.trend_window_sec) || ~isfinite(params.trend_window_sec) || params.trend_window_sec < 0
    error('params.trend_window_sec must be a nonnegative scalar.');
end
end


function [E, info] = local_load_event_density_input(input_data)
info = struct();
info.source_file = '';

if ischar(input_data) || isstring(input_data)
    input_file = char(input_data);
    if exist(input_file, 'file') ~= 2
        error('Input file not found: %s', input_file);
    end
    info.source_file = input_file;
    S = load(input_file, 'E');
    if ~isfield(S, 'E')
        error('Saved file does not contain variable E: %s', input_file);
    end
    E = S.E;
else
    E = input_data;
end

if ~isstruct(E)
    error('input_data must resolve to an event-density struct.');
end
end


function [density_to_plot, raw_density] = local_resolve_density_series(E, params)
raw_density = [];

if params.use_smoothed && isfield(E, 'smoothed_density_mean') && ~isempty(E.smoothed_density_mean)
    density_to_plot = double(E.smoothed_density_mean);
else
    density_to_plot = [];
end

if isempty(density_to_plot)
    if ~isfield(E, 'density_mean') || isempty(E.density_mean)
        error('Event-density struct does not contain density_mean or smoothed_density_mean.');
    end
    density_to_plot = double(E.density_mean);
end

if isfield(E, 'density_mean') && ~isempty(E.density_mean)
    raw_density = double(E.density_mean);
end

if size(density_to_plot, 1) == 1
    density_to_plot = density_to_plot(:);
end
if ~isempty(raw_density) && size(raw_density, 1) == 1
    raw_density = raw_density(:);
end

if ~isempty(raw_density) && size(raw_density, 1) ~= size(density_to_plot, 1)
    raw_density = [];
end
end


function density_to_plot = local_apply_trend_smoothing(density_to_plot, E, t_centers, params)
if isempty(density_to_plot) || params.trend_window_sec <= 0
    return;
end

if isfield(E, 'bin_sec') && ~isempty(E.bin_sec)
    bin_sec = double(E.bin_sec);
elseif numel(t_centers) > 1
    dt = diff(double(t_centers(:)));
    dt = dt(isfinite(dt) & dt > 0);
    if isempty(dt)
        return;
    end
    bin_sec = median(dt);
else
    return;
end

window_bins = max(1, round(double(params.trend_window_sec) / bin_sec));
if window_bins <= 1
    return;
end

density_to_plot = movmean(density_to_plot, window_bins, 1, 'Endpoints', 'shrink');
end


function t_centers = local_resolve_time_axis(E, n_bins)
if isfield(E, 't_centers') && ~isempty(E.t_centers)
    t_centers = double(E.t_centers(:));
elseif isfield(E, 'bin_sec') && ~isempty(E.bin_sec)
    t_centers = ((0:n_bins-1)' + 0.5) * double(E.bin_sec);
else
    error('Event-density struct does not contain t_centers or bin_sec.');
end

if numel(t_centers) ~= n_bins
    error('t_centers length does not match density series length.');
end
end


function [t_plot, unit_label, time_scale] = local_convert_time_axis(t_sec, time_unit)
total_duration_sec = t_sec(end) - t_sec(1);
time_unit = lower(char(string(time_unit)));

switch time_unit
    case 'auto'
        if total_duration_sec >= 2 * 3600
            time_unit = 'hours';
        elseif total_duration_sec >= 10 * 60
            time_unit = 'minutes';
        else
            time_unit = 'seconds';
        end
    case {'sec', 'second', 'seconds'}
        time_unit = 'seconds';
    case {'min', 'minute', 'minutes'}
        time_unit = 'minutes';
    case {'hr', 'hour', 'hours'}
        time_unit = 'hours';
    otherwise
        error('Unsupported time_unit: %s', time_unit);
end

switch time_unit
    case 'seconds'
        time_scale = 1;
        unit_label = 's';
    case 'minutes'
        time_scale = 60;
        unit_label = 'min';
    case 'hours'
        time_scale = 3600;
        unit_label = 'h';
end

t_plot = t_sec(:) / time_scale;
end


function [boundary_t, session_mid_t, session_labels] = local_resolve_session_time_info(E, time_scale)
boundary_t = zeros(0, 1);
session_mid_t = zeros(0, 1);
session_labels = {};

if ~isfield(E, 'session_lengths') || ~isfield(E, 'session_dx') || ...
        isempty(E.session_lengths) || isempty(E.session_dx)
    return;
end

session_lengths = double(E.session_lengths(:));
session_dx = double(E.session_dx(:));
if numel(session_lengths) ~= numel(session_dx)
    return;
end

[~, session_start_time, session_end_time] = io_utils.build_global_time_axis_from_sessions( ...
    session_lengths, session_dx);
if numel(session_start_time) > 1
    boundary_t = session_start_time(2:end) / time_scale;
end
session_mid_t = ((session_start_time + session_end_time) / 2) / time_scale;

if isfield(E, 'session_ids') && ~isempty(E.session_ids) && numel(E.session_ids) == numel(session_mid_t)
    session_ids = double(E.session_ids(:));
else
    session_ids = (1:numel(session_mid_t)).';
end
session_labels = arrayfun(@(x) sprintf('S%d', x), session_ids, 'UniformOutput', false);
end


function band_labels = local_resolve_band_labels(E, n_bands)
if isfield(E, 'band_labels') && ~isempty(E.band_labels)
    band_labels = cellstr(string(E.band_labels(:)));
else
    band_labels = arrayfun(@(i) sprintf('band%d', i), 1:n_bands, 'UniformOutput', false);
end

if numel(band_labels) < n_bands
    for i = numel(band_labels)+1:n_bands
        band_labels{i} = sprintf('band%d', i); %#ok<AGROW>
    end
end
band_labels = band_labels(:).';
end


function dataset_label = local_resolve_dataset_label(E, info)
if isfield(E, 'dataset_id') && ~isempty(E.dataset_id)
    dataset_label = char(string(E.dataset_id));
elseif isfield(E, 'file_stem') && ~isempty(E.file_stem)
    dataset_label = char(string(E.file_stem));
elseif isfield(E, 'source_event_file') && ~isempty(E.source_event_file)
    [~, dataset_label] = fileparts(char(string(E.source_event_file)));
elseif isfield(info, 'source_file') && ~isempty(info.source_file)
    [~, dataset_label] = fileparts(info.source_file);
else
    dataset_label = 'dataset';
end
end


function title_text = local_build_title(dataset_label, E, total_duration_sec, params)
duration_hr = double(total_duration_sec) / 3600;
if isfield(E, 'bin_sec') && ~isempty(E.bin_sec)
    bin_text = sprintf('bin %.2gs', double(E.bin_sec));
else
    bin_text = 'bin n/a';
end
if isfield(E, 'smooth_sigma_sec') && ~isempty(E.smooth_sigma_sec)
    smooth_text = sprintf('smooth %.2gs', double(E.smooth_sigma_sec));
else
    smooth_text = 'smooth n/a';
end

title_text = sprintf('%s pipeline2 event density long-term trend | %.2f h | %s | %s', ...
    dataset_label, duration_hr, bin_text, smooth_text);
if params.trend_window_sec > 0
    title_text = sprintf('%s | trend movmean %.0fs', title_text, double(params.trend_window_sec));
end
if ~isempty(params.title_suffix)
    title_text = sprintf('%s | %s', title_text, char(string(params.title_suffix)));
end
end


function total_duration_sec = local_get_total_duration_sec(E, t_centers)
if isfield(E, 't_edges') && ~isempty(E.t_edges)
    edges = double(E.t_edges(:));
    total_duration_sec = edges(end) - edges(1);
elseif isfield(E, 'session_lengths') && isfield(E, 'session_dx') && ...
        ~isempty(E.session_lengths) && ~isempty(E.session_dx)
    total_duration_sec = sum(double(E.session_lengths(:)) .* double(E.session_dx(:)));
else
    total_duration_sec = t_centers(end) - t_centers(1);
    if isfield(E, 'bin_sec') && ~isempty(E.bin_sec)
        total_duration_sec = total_duration_sec + double(E.bin_sec);
    end
end
end


function figure_prefix = local_build_figure_prefix(E, info, params)
if isfield(params, 'figure_prefix') && ~isempty(params.figure_prefix)
    figure_prefix = char(string(params.figure_prefix));
    return;
end

if isfield(E, 'file_stem') && ~isempty(E.file_stem)
    base = char(string(E.file_stem));
elseif isfield(E, 'dataset_id') && ~isempty(E.dataset_id)
    base = char(string(E.dataset_id));
elseif isfield(info, 'source_file') && ~isempty(info.source_file)
    [~, base] = fileparts(info.source_file);
else
    base = 'dataset';
end

figure_prefix = sprintf('%s_pipeline2_event_density_long_term_trend', base);
if params.trend_window_sec > 0
    figure_prefix = sprintf('%s_movmean%gs', figure_prefix, double(params.trend_window_sec));
end
end


function output_files = local_save_figure(hfig, output_base, save_formats, resolution)
save_formats = cellstr(string(save_formats(:)));
output_files = strings(0, 1);

for i = 1:numel(save_formats)
    fmt = lower(save_formats{i});
    switch fmt
        case 'png'
            output_file = [output_base, '.png'];
            exportgraphics(hfig, output_file, 'Resolution', resolution);
        case 'pdf'
            output_file = [output_base, '.pdf'];
            exportgraphics(hfig, output_file, 'Resolution', resolution, 'ContentType', 'vector');
        case 'fig'
            output_file = [output_base, '.fig'];
            savefig(hfig, output_file);
        otherwise
            error('Unsupported save format: %s', fmt);
    end
    output_files(end+1, 1) = string(output_file); %#ok<AGROW>
end
end


function c = local_lighten_color(c0, mix_with_white)
mix_with_white = max(0, min(1, double(mix_with_white)));
c = double(c0) * (1 - mix_with_white) + mix_with_white;
end


function value = local_get_optional_string_field(S, field_name)
if isfield(S, field_name) && ~isempty(S.(field_name))
    value = string(S.(field_name));
else
    value = "";
end
end
