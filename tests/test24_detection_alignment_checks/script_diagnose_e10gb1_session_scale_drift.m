% Diagnose per-session scale drift for the current E10gb1 raw-data root.
%
% Outputs:
%   tests/test24_detection_alignment_checks/outputs_e10gb1_session_scale_drift/*.csv
%   tests/test24_detection_alignment_checks/outputs_e10gb1_session_scale_drift/*.png

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_script_dir));
addpath(genpath(repo_root));

close all force;

cfg = cfg_E10gb1();
output_dir = fullfile(this_script_dir, 'outputs_e10gb1_session_scale_drift');
if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
end

[session_ids, session_group_labels] = local_collect_included_sessions(cfg);
n_sessions = numel(session_ids);
selected_channels = cfg.channels.selected_all(:);
n_channels = numel(selected_channels);
channel_labels = strings(n_channels, 1);
for c = 1:n_channels
    channel_labels(c) = sprintf('%s_%d', cfg.channels.sites{selected_channels(c)}, selected_channels(c));
end

data_dir = fullfile(cfg.raw_data_root, cfg.data_subfolder);
mean_mat = zeros(n_sessions, n_channels);
std_mat = zeros(n_sessions, n_channels);
length_vec = zeros(n_sessions, 1);
dx_vec = zeros(n_sessions, 1);

for s = 1:n_sessions
    sid = session_ids(s);
    file_path = fullfile(data_dir, sprintf('%s_%04d_blp.mat', cfg.file_stem, sid));
    S = load(file_path, 'blp');
    x = double(S.blp.dat(:, selected_channels, 1));

    mean_mat(s, :) = mean(x, 1);
    std_mat(s, :) = std(x, 0, 1);
    length_vec(s) = size(x, 1);
    dx_vec(s) = double(S.blp.dx);
end

session_stats_table = local_build_session_stats_table( ...
    session_ids, session_group_labels, length_vec, dx_vec, channel_labels, mean_mat, std_mat);
trend_summary_table = local_build_trend_summary_table( ...
    session_ids, session_group_labels, channel_labels, mean_mat, std_mat);

session_plot_file = fullfile(output_dir, 'e10gb1_session_std_by_channel.png');
summary_plot_file = fullfile(output_dir, 'e10gb1_session_std_summary.png');

local_plot_channel_std(session_ids, session_group_labels, channel_labels, std_mat, session_plot_file);
local_plot_summary_std(session_ids, session_group_labels, channel_labels, std_mat, summary_plot_file);

save_file = fullfile(output_dir, 'e10gb1_session_scale_drift.mat');
save(save_file, ...
    'cfg', ...
    'session_ids', ...
    'session_group_labels', ...
    'channel_labels', ...
    'length_vec', ...
    'dx_vec', ...
    'mean_mat', ...
    'std_mat', ...
    'session_stats_table', ...
    'trend_summary_table', ...
    'session_plot_file', ...
    'summary_plot_file', ...
    '-v7.3');

writetable(session_stats_table, fullfile(output_dir, 'session_stats_table.csv'));
writetable(trend_summary_table, fullfile(output_dir, 'trend_summary_table.csv'));

disp('Per-channel trend summary:');
disp(trend_summary_table);

fprintf('Saved session-scale drift outputs to:\n  %s\n', output_dir);


function [session_ids, session_group_labels] = local_collect_included_sessions(cfg)
session_ids = zeros(0, 1);
session_group_labels = strings(0, 1);

for i = 1:numel(cfg.sessions)
    if ~cfg.sessions(i).include
        continue;
    end

    ids = double(cfg.sessions(i).session_id(:));
    labels = repmat(string(cfg.sessions(i).notes), numel(ids), 1);

    session_ids = [session_ids; ids]; %#ok<AGROW>
    session_group_labels = [session_group_labels; labels]; %#ok<AGROW>
end
end


function tbl = local_build_session_stats_table(session_ids, session_group_labels, length_vec, dx_vec, channel_labels, mean_mat, std_mat)
n_sessions = numel(session_ids);
n_channels = numel(channel_labels);
rows = cell(n_sessions * n_channels, 1);
idx = 0;

for s = 1:n_sessions
    for c = 1:n_channels
        idx = idx + 1;
        rows{idx} = table( ...
            s, ...
            session_ids(s), ...
            session_group_labels(s), ...
            length_vec(s), ...
            dx_vec(s), ...
            channel_labels(c), ...
            mean_mat(s, c), ...
            std_mat(s, c), ...
            'VariableNames', { ...
                'session_order', ...
                'session_id', ...
                'session_group', ...
                'n_samples', ...
                'dx', ...
                'channel_label', ...
                'channel_mean', ...
                'channel_std'});
    end
end

tbl = vertcat(rows{:});
end


function tbl = local_build_trend_summary_table(session_ids, session_group_labels, channel_labels, mean_mat, std_mat)
n_channels = numel(channel_labels);
rows = cell(n_channels, 1);

session_order = transpose(1:numel(session_ids));
is_polar = session_group_labels == "polar";
is_spont = session_group_labels == "spont";

for c = 1:n_channels
    mean_vals = mean_mat(:, c);
    std_vals = std_mat(:, c);

    rho_all = local_safe_corr(session_order, std_vals);
    rho_spont = local_safe_corr(session_order(is_spont), std_vals(is_spont));

    rows{c} = table( ...
        channel_labels(c), ...
        min(mean_vals), ...
        max(mean_vals), ...
        max(mean_vals) - min(mean_vals), ...
        max(abs(mean_vals)), ...
        min(std_vals), ...
        max(std_vals), ...
        max(std_vals) / min(std_vals), ...
        mean(std_vals(is_polar)), ...
        mean(std_vals(is_spont)), ...
        mean(std_vals(is_spont)) / mean(std_vals(is_polar)), ...
        rho_all, ...
        rho_spont, ...
        'VariableNames', { ...
            'channel_label', ...
            'mean_min', ...
            'mean_max', ...
            'mean_range', ...
            'mean_absmax', ...
            'std_min', ...
            'std_max', ...
            'std_ratio_max_over_min', ...
            'polar_mean_std', ...
            'spont_mean_std', ...
            'spont_to_polar_std_ratio', ...
            'std_session_order_corr_all', ...
            'std_session_order_corr_spont_only'});
end

tbl = vertcat(rows{:});
end


function rho = local_safe_corr(x, y)
rho = NaN;
if numel(x) < 2 || numel(y) < 2
    return;
end
if all(abs(y - y(1)) < eps)
    return;
end
rho = corr(x(:), y(:), 'type', 'Spearman', 'rows', 'complete');
end


function local_plot_channel_std(session_ids, session_group_labels, channel_labels, std_mat, save_file)
n_channels = numel(channel_labels);
fig = figure('Color', 'w', 'Position', [100, 100, 1500, 900]);
tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

polar_color = [0.00, 0.45, 0.74];
spont_color = [0.85, 0.33, 0.10];

for c = 1:n_channels
    nexttile;
    hold on;

    y = std_mat(:, c);
    plot(session_ids, y, '-', 'Color', [0.65, 0.65, 0.65], 'LineWidth', 1);

    polar_mask = session_group_labels == "polar";
    spont_mask = session_group_labels == "spont";

    scatter(session_ids(polar_mask), y(polar_mask), 36, polar_color, 'filled');
    scatter(session_ids(spont_mask), y(spont_mask), 36, spont_color, 'filled');

    title(strrep(channel_labels(c), '_', '\_'), 'Interpreter', 'tex');
    xlabel('Session ID');
    ylabel('Std');
    grid on;
    box off;
end

lgd = legend({'session order', 'polar', 'spont'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south';
sgtitle('E10gb1 per-session raw-channel std by selected channel');

exportgraphics(fig, save_file, 'Resolution', 200);
close(fig);
end


function local_plot_summary_std(session_ids, session_group_labels, channel_labels, std_mat, save_file)
fig = figure('Color', 'w', 'Position', [120, 120, 1200, 700]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

polar_color = [0.00, 0.45, 0.74];
spont_color = [0.85, 0.33, 0.10];
hp_mask = startsWith(channel_labels, "hp_");
pl_mask = startsWith(channel_labels, "pl_");

session_polar = session_group_labels == "polar";
session_spont = session_group_labels == "spont";

nexttile;
hold on;
hp_mean_std = mean(std_mat(:, hp_mask), 2);
plot(session_ids, hp_mean_std, '-', 'Color', [0.65, 0.65, 0.65], 'LineWidth', 1);
scatter(session_ids(session_polar), hp_mean_std(session_polar), 42, polar_color, 'filled');
scatter(session_ids(session_spont), hp_mean_std(session_spont), 42, spont_color, 'filled');
title('Mean session std across hp channels');
xlabel('Session ID');
ylabel('Mean std');
grid on;
box off;

nexttile;
hold on;
pl_mean_std = mean(std_mat(:, pl_mask), 2);
plot(session_ids, pl_mean_std, '-', 'Color', [0.65, 0.65, 0.65], 'LineWidth', 1);
scatter(session_ids(session_polar), pl_mean_std(session_polar), 42, polar_color, 'filled');
scatter(session_ids(session_spont), pl_mean_std(session_spont), 42, spont_color, 'filled');
title('Mean session std across pl channels');
xlabel('Session ID');
ylabel('Mean std');
grid on;
box off;

lgd = legend({'session order', 'polar', 'spont'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south';
sgtitle('E10gb1 per-session raw-channel std summary');

exportgraphics(fig, save_file, 'Resolution', 200);
close(fig);
end
