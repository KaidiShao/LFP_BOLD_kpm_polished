% Run lagged spike-vs-LFP-residual correlation for current P6 outputs.
%
% Optional variables before run(...):
%   dataset_stems      = {'e10gb1','e10fV1','e10gh1','f12m01','e10gw1'};
%   max_lag_bins       = 80;  % about +/-2 s when spike dx is 0.025 s
%   lag_features       = {'abs_rms', 'abs_mean'};
%   force_recompute    = false;
%   continue_on_error  = true;
%   dry_run            = false;
%   max_tasks          = [];

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));
set(groot, 'defaultFigureVisible', 'off');

if ~exist('dataset_stems', 'var') || isempty(dataset_stems)
    dataset_stems = {'e10gb1', 'e10fV1', 'e10gh1', 'f12m01', 'e10gw1'};
end
if ~exist('max_lag_bins', 'var') || isempty(max_lag_bins)
    max_lag_bins = 80;
end
if ~exist('lag_features', 'var') || isempty(lag_features)
    lag_features = {'abs_rms', 'abs_mean'};
end
if ~exist('force_recompute', 'var') || isempty(force_recompute)
    force_recompute = false;
end
if ~exist('continue_on_error', 'var') || isempty(continue_on_error)
    continue_on_error = true;
end
if ~exist('dry_run', 'var') || isempty(dry_run)
    dry_run = false;
end
if ~exist('max_tasks', 'var')
    max_tasks = [];
end

dataset_stems = cellstr(string(dataset_stems(:)).');
lag_features = cellstr(string(lag_features(:)).');
max_lag_bins = double(max_lag_bins);
force_recompute = logical(force_recompute);
continue_on_error = logical(continue_on_error);
dry_run = logical(dry_run);
if ~isempty(max_tasks)
    max_tasks = double(max_tasks);
end

processed_root = io_project.get_project_processed_root();
source_stage = io_project.get_pipeline_stage_name(6, 'spkt_residual_cross_correlation');
target_stage = io_project.get_pipeline_stage_name(6, 'spkt_residual_lagged_cross_correlation');
figure_stage = io_project.get_pipeline_stage_name(6, 'figures_spkt_residual_lagged_cross_correlation');

fprintf('Processed root: %s\n', processed_root);
fprintf('Source stage : %s\n', source_stage);
fprintf('Target stage : %s\n', target_stage);
fprintf('Figure stage : %s\n', figure_stage);
fprintf('Datasets     : %s\n', strjoin(dataset_stems, ', '));
fprintf('Features     : %s\n', strjoin(lag_features, ', '));
fprintf('Max lag bins : %d\n', max_lag_bins);
fprintf('Dry run      : %d\n\n', dry_run);

tasks = {};
for i_ds = 1:numel(dataset_stems)
    dataset_stem = dataset_stems{i_ds};
    source_dir = fullfile(processed_root, dataset_stem, source_stage);
    if exist(source_dir, 'dir') ~= 7
        fprintf('[skip] missing source dir: %s\n', source_dir);
        continue;
    end

    L = dir(fullfile(source_dir, '*.mat'));
    L = L(~contains({L.name}, '_lagged', 'IgnoreCase', true));
    [~, order] = sort({L.name});
    L = L(order);

    save_dir = fullfile(processed_root, dataset_stem, target_stage);
    figure_dir = fullfile(processed_root, dataset_stem, figure_stage);
    for i_file = 1:numel(L)
        source_file = fullfile(L(i_file).folder, L(i_file).name);
        [~, source_stem] = fileparts(source_file);
        save_tag = sprintf('%s_lagged_pm%dbins', source_stem, max_lag_bins);
        out_mat = fullfile(save_dir, [local_filename_safe(save_tag), '.mat']);
        out_png = fullfile(figure_dir, [local_filename_safe(save_tag), '_overview.png']);

        tasks{end + 1, 1} = struct( ... %#ok<SAGROW>
            'dataset_stem', dataset_stem, ...
            'source_file', source_file, ...
            'save_dir', save_dir, ...
            'figure_dir', figure_dir, ...
            'save_tag', save_tag, ...
            'out_mat', out_mat, ...
            'out_png', out_png);
    end
end

fprintf('Discovered %d lagged SPKT task(s).\n\n', numel(tasks));
if ~isempty(max_tasks) && isfinite(max_tasks) && max_tasks > 0
    tasks = tasks(1:min(numel(tasks), round(max_tasks)));
    fprintf('Limiting to first %d task(s).\n\n', numel(tasks));
end

if dry_run
    for i_task = 1:numel(tasks)
        task = tasks{i_task};
        fprintf('[dry-run %02d/%02d] %s\n', i_task, numel(tasks), task.dataset_stem);
        fprintf('  source: %s\n', task.source_file);
        fprintf('  out   : %s\n', task.out_mat);
        fprintf('  fig   : %s\n', task.out_png);
    end
    return;
end

ok_count = 0;
skip_count = 0;
fail_count = 0;

for i_task = 1:numel(tasks)
    task = tasks{i_task};
    fprintf('[%02d/%02d] %s\n', i_task, numel(tasks), task.source_file);

    if ~force_recompute && exist(task.out_mat, 'file') == 2 && exist(task.out_png, 'file') == 2
        fprintf('  skip existing:\n    %s\n    %s\n\n', task.out_mat, task.out_png);
        skip_count = skip_count + 1;
        continue;
    end

    try
        if exist(task.save_dir, 'dir') ~= 7
            mkdir(task.save_dir);
        end
        if exist(task.figure_dir, 'dir') ~= 7
            mkdir(task.figure_dir);
        end

        params = struct();
        params.features = lag_features;
        params.max_lag_bins = max_lag_bins;
        params.top_n_rows = 200;
        params.save_results = true;
        params.save_dir = task.save_dir;
        params.save_tag = task.save_tag;
        params.verbose = true;

        out = compute_spike_residual_lagged_correlation(task.source_file, params);
        local_export_lagged_overview(out, task.figure_dir, local_filename_safe(task.save_tag));

        ok_count = ok_count + 1;
        fprintf('  OK -> %s\n\n', out.save_paths.main_mat);
    catch ME
        fail_count = fail_count + 1;
        fprintf(2, '  FAILED: %s\n', ME.message);
        if ~continue_on_error
            rethrow(ME);
        end
        fprintf('\n');
    end
end

fprintf('Done. OK: %d | skipped: %d | failed: %d\n', ok_count, skip_count, fail_count);


function local_export_lagged_overview(out, figure_dir, save_tag)
if exist(figure_dir, 'dir') ~= 7
    mkdir(figure_dir);
end

feature_primary = out.features{1};
if numel(out.features) >= 2
    feature_secondary = out.features{2};
else
    feature_secondary = out.features{1};
end

channel_labels = cellstr(string(out.channel_table.channel_label));
mode_labels = compose('m%d', double(out.mode_table.mode_rank));
n_channels = height(out.channel_table);
n_modes = height(out.mode_table);

peak_abs_primary = out.lag_corr.([feature_primary, '_peak']).peak_abs_corr;
peak_lag_primary = out.lag_corr.([feature_primary, '_peak']).peak_lag_sec;
peak_abs_secondary = out.lag_corr.([feature_secondary, '_peak']).peak_abs_corr;
peak_lag_secondary = out.lag_corr.([feature_secondary, '_peak']).peak_lag_sec;

mag_vals = [peak_abs_primary(:); peak_abs_secondary(:)];
mag_vals = mag_vals(isfinite(mag_vals));
if isempty(mag_vals)
    mag_clim = [0, 1];
else
    mag_clim = [0, max(mag_vals)];
    if mag_clim(2) == 0
        mag_clim(2) = 1;
    end
end

lag_vals = [peak_lag_primary(:); peak_lag_secondary(:)];
lag_vals = lag_vals(isfinite(lag_vals));
if isempty(lag_vals)
    lag_clim = [-max(abs(out.lag_sec)), max(abs(out.lag_sec))];
else
    lag_lim = max(abs(lag_vals));
    if lag_lim == 0
        lag_lim = max(abs(out.lag_sec));
    end
    lag_clim = [-lag_lim, lag_lim];
end

hfig = figure('Color', 'w', 'Position', [80, 60, 1850, 980]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(1);
imagesc(peak_abs_primary);
axis tight;
set(gca, 'YTick', 1:n_channels, 'YTickLabel', channel_labels, ...
    'XTick', 1:n_modes, 'XTickLabel', mode_labels, ...
    'TickLabelInterpreter', 'none');
xtickangle(0);
xlabel('Residual Mode Rank');
ylabel('Spike Channel');
title(sprintf('Peak |corr|: %s', strrep(feature_primary, '_', '\_')));
colormap(gca, parula(256));
colorbar;
caxis(mag_clim);

nexttile(2);
imagesc(peak_lag_primary);
axis tight;
set(gca, 'YTick', 1:n_channels, 'YTickLabel', channel_labels, ...
    'XTick', 1:n_modes, 'XTickLabel', mode_labels, ...
    'TickLabelInterpreter', 'none');
xtickangle(0);
xlabel('Residual Mode Rank');
ylabel('Spike Channel');
title(sprintf('Peak lag (s): %s', strrep(feature_primary, '_', '\_')));
colormap(gca, local_blue_white_red(256));
colorbar;
caxis(lag_clim);

nexttile(3, [2 1]);
local_plot_top_curves(out, feature_primary, 6);
title(sprintf('Top lag curves: %s', strrep(feature_primary, '_', '\_')));

nexttile(4);
imagesc(peak_abs_secondary);
axis tight;
set(gca, 'YTick', 1:n_channels, 'YTickLabel', channel_labels, ...
    'XTick', 1:n_modes, 'XTickLabel', mode_labels, ...
    'TickLabelInterpreter', 'none');
xtickangle(0);
xlabel('Residual Mode Rank');
ylabel('Spike Channel');
title(sprintf('Peak |corr|: %s', strrep(feature_secondary, '_', '\_')));
colormap(gca, parula(256));
colorbar;
caxis(mag_clim);

nexttile(5);
imagesc(peak_lag_secondary);
axis tight;
set(gca, 'YTick', 1:n_channels, 'YTickLabel', channel_labels, ...
    'XTick', 1:n_modes, 'XTickLabel', mode_labels, ...
    'TickLabelInterpreter', 'none');
xtickangle(0);
xlabel('Residual Mode Rank');
ylabel('Spike Channel');
title(sprintf('Peak lag (s): %s', strrep(feature_secondary, '_', '\_')));
colormap(gca, local_blue_white_red(256));
colorbar;
caxis(lag_clim);

sgtitle(sprintf(['Lagged spike vs LFP residual correlation\n', ...
    'Positive lag means spike leads residual | max lag = %.2f s'], max(abs(out.lag_sec))), ...
    'FontWeight', 'bold');

png_file = fullfile(figure_dir, [save_tag, '_overview.png']);
fig_file = fullfile(figure_dir, [save_tag, '_overview.fig']);
exportgraphics(hfig, png_file, 'Resolution', 220);
savefig(hfig, fig_file);
if isgraphics(hfig)
    close(hfig);
end

fprintf('  Saved overview PNG: %s\n', png_file);
fprintf('  Saved overview FIG: %s\n', fig_file);
end


function local_plot_top_curves(out, feature_name, top_n)
T = out.peak_table;
mask = strcmp(string(T.feature), string(feature_name)) & isfinite(T.peak_abs_corr);
Tf = sortrows(T(mask, :), 'peak_abs_corr', 'descend');
top_n = min(top_n, height(Tf));

hold on;
grid on;
box on;

if top_n == 0
    text(0.5, 0.5, 'No finite lagged correlations found', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    xlabel('Lag (s)');
    ylabel('Correlation');
    return;
end

curve_cube = out.lag_corr.(feature_name);
colors = lines(top_n);

for i = 1:top_n
    ch_row = find(out.channel_table.channel_index == Tf.channel_index(i), 1, 'first');
    mode_row = find(out.mode_table.mode_rank == Tf.mode_rank(i), 1, 'first');
    if isempty(ch_row) || isempty(mode_row)
        continue;
    end

    curve = squeeze(curve_cube(ch_row, mode_row, :));
    plot(out.lag_sec, curve, 'LineWidth', 1.6, 'Color', colors(i, :), ...
        'DisplayName', sprintf('%s | m%d | peak %.3fs', ...
        char(Tf.channel_label(i)), double(Tf.mode_rank(i)), double(Tf.peak_lag_sec(i))));
    plot(double(Tf.peak_lag_sec(i)), double(Tf.peak_corr(i)), 'o', ...
        'MarkerSize', 7, 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :), ...
        'HandleVisibility', 'off');
end

xline(0, 'k--', 'LineWidth', 1.0);
yline(0, 'k-', 'LineWidth', 0.8);
xlabel('Lag (s)');
ylabel('Correlation');
legend('Location', 'eastoutside', 'Interpreter', 'none');
end


function cmap = local_blue_white_red(n)
if nargin < 1
    n = 256;
end
half_n = ceil(n / 2);
blue = [0.13, 0.35, 0.75];
white = [1, 1, 1];
red = [0.76, 0.16, 0.12];

left = [linspace(blue(1), white(1), half_n).', ...
    linspace(blue(2), white(2), half_n).', ...
    linspace(blue(3), white(3), half_n).'];
right = [linspace(white(1), red(1), n - half_n + 1).', ...
    linspace(white(2), red(2), n - half_n + 1).', ...
    linspace(white(3), red(3), n - half_n + 1).'];
cmap = [left; right(2:end, :)];
end


function out = local_filename_safe(name_in)
out = regexprep(char(name_in), '[^a-zA-Z0-9_\-]+', '_');
out = regexprep(out, '_+', '_');
out = regexprep(out, '^_+|_+$', '');
if isempty(out)
    out = 'untitled';
end
end
