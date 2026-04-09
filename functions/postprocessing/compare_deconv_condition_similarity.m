function [fig, sim_info] = compare_deconv_condition_similarity(deconv_results, cfg)
%COMPARE_DECONV_CONDITION_SIMILARITY Compare heatmap values across deconv conditions.
%
%   [fig, sim_info] = compare_deconv_condition_similarity(deconv_results, cfg)
%
% Inputs
%   deconv_results  nested struct produced by script_visualize_edmd_outputs_windowed.m
%   cfg             optional struct
%     .method_order         cellstr order for method dimension
%     .lambda_source_order  cellstr order for lambda-source dimension
%     .scope_order          cellstr order for normalisation dimension
%     .order_mode           'scope_first' (default) or 'method_first'
%     .compare_fields       cellstr of plot_view fields (default {'abs_all','real_all'})
%     .title_prefix         figure title prefix (default 'Deconv condition similarity')
%
% Output
%   fig       figure handle
%   sim_info  struct containing labels, matrices, and pairwise summaries

if nargin < 2 || isempty(cfg)
    cfg = struct();
end

cfg = local_set_default(cfg, 'method_order', fieldnames(deconv_results));
cfg = local_set_default(cfg, 'lambda_source_order', {});
cfg = local_set_default(cfg, 'scope_order', {});
cfg = local_set_default(cfg, 'order_mode', 'scope_first');
cfg = local_set_default(cfg, 'compare_fields', {'abs_all', 'real_all'});
cfg = local_set_default(cfg, 'title_prefix', 'Deconv condition similarity');

cfg.method_order = local_make_valid_name_list(cfg.method_order);
cfg.lambda_source_order = local_make_valid_name_list(cfg.lambda_source_order);
cfg.scope_order = local_make_valid_name_list(cfg.scope_order);

[entries, labels] = local_collect_entries(deconv_results, cfg);
n_cond = numel(entries);

if n_cond < 2
    error('Need at least two deconv conditions to compare.');
end

n_fields = numel(cfg.compare_fields);
sim_info = struct();
sim_info.cfg = cfg;
sim_info.labels = labels;
sim_info.entries = entries;
sim_info.n_conditions = n_cond;

fig = figure('Color', 'w', 'Name', cfg.title_prefix);
tl = tiledlayout(fig, 2, n_fields, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, cfg.title_prefix, 'Interpreter', 'none');

for i_field = 1:n_fields
    field_name = cfg.compare_fields{i_field};
    [corr_mat, mad_mat] = local_compare_field(entries, field_name);
    sim_info.(field_name) = struct();
    sim_info.(field_name).corr = corr_mat;
    sim_info.(field_name).mad = mad_mat;
    sim_info.(field_name).summary = local_pair_summary(corr_mat, mad_mat, labels);

    nexttile(tl, i_field);
    imagesc(corr_mat);
    axis image;
    set(gca, 'YDir', 'normal', 'XTick', 1:n_cond, 'YTick', 1:n_cond, ...
        'XTickLabel', labels, 'YTickLabel', labels, 'XTickLabelRotation', 45);
    clim([0 1]);
    colormap(gca, parula(256));
    colorbar;
    title(sprintf('%s: corr', local_field_title(field_name)), 'Interpreter', 'none');

    nexttile(tl, n_fields + i_field);
    imagesc(mad_mat);
    axis image;
    set(gca, 'YDir', 'normal', 'XTick', 1:n_cond, 'YTick', 1:n_cond, ...
        'XTickLabel', labels, 'YTickLabel', labels, 'XTickLabelRotation', 45);
    colormap(gca, flipud(gray(256)));
    colorbar;
    title(sprintf('%s: mean |diff|', local_field_title(field_name)), 'Interpreter', 'none');
end

fig.UserData.sim_info = sim_info;
end

function cfg = local_set_default(cfg, name, value)
if ~isfield(cfg, name) || isempty(cfg.(name))
    cfg.(name) = value;
end
end

function [entries, labels] = local_collect_entries(deconv_results, cfg)
entries = {};
labels = {};

method_fields = cfg.method_order;
lambda_fields_default = cfg.lambda_source_order;
scope_fields_default = cfg.scope_order;

switch lower(cfg.order_mode)
    case 'scope_first'
        for i_scope = 1:numel(scope_fields_default)
            scope_field = scope_fields_default{i_scope};
            for i_method = 1:numel(method_fields)
                method_field = method_fields{i_method};
                if ~isfield(deconv_results, method_field)
                    continue;
                end

                lambda_fields = lambda_fields_default;
                if isempty(lambda_fields)
                    lambda_fields = fieldnames(deconv_results.(method_field));
                end

                for i_lambda = 1:numel(lambda_fields)
                    lambda_field = lambda_fields{i_lambda};
                    if ~isfield(deconv_results.(method_field), lambda_field)
                        continue;
                    end
                    if ~isfield(deconv_results.(method_field).(lambda_field), scope_field)
                        continue;
                    end

                    deconv = deconv_results.(method_field).(lambda_field).(scope_field);
                    entries{end+1,1} = deconv; %#ok<AGROW>
                    labels{end+1,1} = sprintf('%s | %s | %s', ...
                        local_label_from_field(method_field), ...
                        local_label_from_field(lambda_field), ...
                        local_label_from_field(scope_field)); %#ok<AGROW>
                end
            end
        end

    case 'method_first'
        for i_method = 1:numel(method_fields)
            method_field = method_fields{i_method};
            if ~isfield(deconv_results, method_field)
                continue;
            end

            lambda_fields = lambda_fields_default;
            if isempty(lambda_fields)
                lambda_fields = fieldnames(deconv_results.(method_field));
            end

            for i_lambda = 1:numel(lambda_fields)
                lambda_field = lambda_fields{i_lambda};
                if ~isfield(deconv_results.(method_field), lambda_field)
                    continue;
                end

                scope_fields = scope_fields_default;
                if isempty(scope_fields)
                    scope_fields = fieldnames(deconv_results.(method_field).(lambda_field));
                end

                for i_scope = 1:numel(scope_fields)
                    scope_field = scope_fields{i_scope};
                    if ~isfield(deconv_results.(method_field).(lambda_field), scope_field)
                        continue;
                    end

                    deconv = deconv_results.(method_field).(lambda_field).(scope_field);
                    entries{end+1,1} = deconv; %#ok<AGROW>
                    labels{end+1,1} = sprintf('%s | %s | %s', ...
                        local_label_from_field(method_field), ...
                        local_label_from_field(lambda_field), ...
                        local_label_from_field(scope_field)); %#ok<AGROW>
                end
            end
        end

    otherwise
        error('Unknown cfg.order_mode = %s. Use ''scope_first'' or ''method_first''.', cfg.order_mode);
end
end

function [corr_mat, mad_mat] = local_compare_field(entries, field_name)
n = numel(entries);
corr_mat = eye(n);
mad_mat = zeros(n);
vecs = cell(n, 1);

for i = 1:n
    if ~isfield(entries{i}, 'plot_view') || ~isfield(entries{i}.plot_view, field_name)
        error('Field plot_view.%s not found in one of the deconv entries.', field_name);
    end
    vecs{i} = entries{i}.plot_view.(field_name)(:);
end

ref_size = size(vecs{1});
for i = 2:n
    if ~isequal(size(vecs{i}), ref_size)
        error('Heatmap field %s has inconsistent sizes across conditions.', field_name);
    end
end

for i = 1:n
    for j = i+1:n
        corr_ij = local_safe_corr(vecs{i}, vecs{j});
        mad_ij = mean(abs(vecs{i} - vecs{j}), 'omitnan');
        corr_mat(i, j) = corr_ij;
        corr_mat(j, i) = corr_ij;
        mad_mat(i, j) = mad_ij;
        mad_mat(j, i) = mad_ij;
    end
end
end

function r = local_safe_corr(x, y)
x = x(:);
y = y(:);
keep = isfinite(x) & isfinite(y);
x = x(keep);
y = y(keep);

if isempty(x) || isempty(y)
    r = nan;
    return;
end

x = x - mean(x);
y = y - mean(y);
den = sqrt(sum(x.^2) * sum(y.^2));
if den <= 0
    r = 1;
else
    r = sum(x .* y) / den;
end
end

function summary = local_pair_summary(corr_mat, mad_mat, labels)
n = size(corr_mat, 1);
mask = triu(true(n), 1);

corr_vals = corr_mat(mask);
mad_vals = mad_mat(mask);

[max_corr, idx_max_corr] = max(corr_vals);
[min_corr, idx_min_corr] = min(corr_vals);
[~, idx_min_mad] = min(mad_vals);
[~, idx_max_mad] = max(mad_vals);

[row_idx, col_idx] = find(mask);

summary = struct();
summary.mean_corr = mean(corr_vals, 'omitnan');
summary.median_corr = median(corr_vals, 'omitnan');
summary.mean_mad = mean(mad_vals, 'omitnan');
summary.median_mad = median(mad_vals, 'omitnan');
summary.most_similar_by_corr = local_pair_record(row_idx(idx_max_corr), col_idx(idx_max_corr), ...
    labels, max_corr, mad_mat);
summary.least_similar_by_corr = local_pair_record(row_idx(idx_min_corr), col_idx(idx_min_corr), ...
    labels, min_corr, mad_mat);
summary.most_similar_by_mad = local_pair_record(row_idx(idx_min_mad), col_idx(idx_min_mad), ...
    labels, corr_mat(row_idx(idx_min_mad), col_idx(idx_min_mad)), mad_mat);
summary.least_similar_by_mad = local_pair_record(row_idx(idx_max_mad), col_idx(idx_max_mad), ...
    labels, corr_mat(row_idx(idx_max_mad), col_idx(idx_max_mad)), mad_mat);
end

function rec = local_pair_record(i, j, labels, corr_val, mad_mat)
rec = struct();
rec.idx = [i, j];
rec.label_i = labels{i};
rec.label_j = labels{j};
rec.corr = corr_val;
rec.mad = mad_mat(i, j);
end

function title_str = local_field_title(field_name)
switch lower(field_name)
    case 'abs_all'
        title_str = 'ALL |u|';
    case 'real_all'
        title_str = 'ALL Re(u)';
    case 'abs_sel'
        title_str = 'SELECTED |u|';
    case 'real_sel'
        title_str = 'SELECTED Re(u)';
    otherwise
        title_str = strrep(field_name, '_', ' ');
end
end

function label = local_label_from_field(field_name)
label = strrep(field_name, '_', ' ');
switch lower(field_name)
    case 'koopman_residual'
        label = 'residual';
    case 'empirical_complex'
        label = 'empirical complex';
    case 'xglobal'
        label = 'global';
end
end

function out = local_make_valid_name_list(names_in)
if isempty(names_in)
    out = names_in;
    return;
end

out = cell(size(names_in));
for i = 1:numel(names_in)
    out{i} = matlab.lang.makeValidName(names_in{i});
end
end
