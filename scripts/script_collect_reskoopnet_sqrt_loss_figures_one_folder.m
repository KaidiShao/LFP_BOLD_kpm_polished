repo_root = fileparts(fileparts(mfilename('fullpath')));
source_root = 'E:\autodl_results';
output_root = fullfile( ...
    repo_root, ...
    'analysis_outputs', ...
    ['reskoopnet_training_sqrt_loss_figures_one_folder_' char(datetime('today', 'Format', 'yyyyMMdd'))] ...
);

include_pdf = true;
include_non_squared_runs = false;

if ~isfolder(source_root)
    error('Source root does not exist: %s', source_root);
end

if ~isfolder(output_root)
    mkdir(output_root);
end

summary_files = dir(fullfile(source_root, '**', '*_summary.mat'));
if isempty(summary_files)
    error('No summary MAT files were found under %s', source_root);
end

records = cell(0, 12);
export_count = 0;
skip_count = 0;

for k = 1:numel(summary_files)
    summary_path = fullfile(summary_files(k).folder, summary_files(k).name);
    status = "skipped";
    reason = "";
    diag_png = "";
    metrics_png = "";
    loss_mode = "";
    requested_loss_mode = "";

    try
        S = load(summary_path, 'EDMD_outputs');
        if ~isfield(S, 'EDMD_outputs')
            error('Variable EDMD_outputs is missing.');
        end
        EDMD_outputs = S.EDMD_outputs;

        loss_mode = local_get_text_field(EDMD_outputs, 'loss_mode', "squared");
        requested_loss_mode = local_get_text_field(EDMD_outputs, 'requested_loss_mode', loss_mode);
        if ~include_non_squared_runs && ~strcmpi(loss_mode, "squared")
            reason = "non_squared_loss_mode";
        else
            dataset_stem = local_get_text_field( ...
                EDMD_outputs, 'dataset_stem', local_guess_dataset_stem(summary_path) ...
            );
            observable_mode = local_get_text_field( ...
                EDMD_outputs, 'observable_mode', local_guess_observable_mode(summary_path) ...
            );
            residual_form = local_get_text_field( ...
                EDMD_outputs, 'residual_form', local_guess_residual_form(summary_path) ...
            );
            run_label = local_get_text_field( ...
                EDMD_outputs, 'run_label', string(summary_files(k).folder) ...
            );
            run_label = local_take_last_path_part(run_label);
            family = local_guess_family(summary_path);

            inner_loss = local_get_numeric_field(EDMD_outputs, 'loss');
            inner_val_loss = local_get_numeric_field(EDMD_outputs, 'val_loss');
            if isempty(inner_loss) && isempty(inner_val_loss)
                reason = "missing_inner_loss";
            else
                outer_epoch = local_get_numeric_field(EDMD_outputs, 'outer_epoch_history');
                outer_train_squared = local_first_numeric_field( ...
                    EDMD_outputs, ...
                    {'outer_train_metric_squared_history', 'outer_train_metric_history'} ...
                );
                outer_val_squared = local_first_numeric_field( ...
                    EDMD_outputs, ...
                    {'outer_val_metric_squared_history', 'outer_val_metric_history'} ...
                );
                outer_cond = local_get_numeric_field(EDMD_outputs, 'outer_eigvec_cond_history');

                inner_train_squared = local_first_numeric_field( ...
                    EDMD_outputs, ...
                    {'inner_train_metric_squared_history', 'loss'} ...
                );
                inner_val_squared = local_first_numeric_field( ...
                    EDMD_outputs, ...
                    {'inner_val_metric_squared_history', 'val_loss'} ...
                );
                inner_train_per_dim = local_get_numeric_field( ...
                    EDMD_outputs, 'inner_train_metric_per_dim_history' ...
                );
                inner_val_per_dim = local_get_numeric_field( ...
                    EDMD_outputs, 'inner_val_metric_per_dim_history' ...
                );

                if isempty(outer_epoch)
                    outer_epoch = 1:max(numel(outer_train_squared), numel(outer_val_squared));
                end

                sqrt_inner_loss = local_sqrt_nonnegative(inner_loss);
                sqrt_inner_val_loss = local_sqrt_nonnegative(inner_val_loss);
                sqrt_outer_train = local_sqrt_nonnegative(outer_train_squared);
                sqrt_outer_val = local_sqrt_nonnegative(outer_val_squared);
                sqrt_inner_squared = local_sqrt_nonnegative(inner_train_squared);
                sqrt_inner_val_squared = local_sqrt_nonnegative(inner_val_squared);
                sqrt_inner_per_dim = local_sqrt_nonnegative(inner_train_per_dim);
                sqrt_inner_val_per_dim = local_sqrt_nonnegative(inner_val_per_dim);

                base_tag = local_sanitize_filename(strjoin([ ...
                    family, dataset_stem, residual_form, observable_mode, run_label ...
                ], '__'));

                fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1200 820]);
                t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
                title(t, sprintf('%s %s %s sqrt squared residual diagnostics', ...
                    dataset_stem, observable_mode, residual_form), 'Interpreter', 'none');

                ax = nexttile(t, 1);
                hold(ax, 'on');
                local_plot_series(ax, 1:numel(sqrt_inner_loss), sqrt_inner_loss, [0 0.447 0.741], 'train');
                local_plot_series(ax, 1:numel(sqrt_inner_val_loss), sqrt_inner_val_loss, [0.85 0.325 0.098], 'val');
                xlabel(ax, 'inner fit epoch');
                ylabel(ax, 'sqrt(loss)');
                title(ax, 'Inner RMS residual');
                grid(ax, 'on');
                legend(ax, 'Location', 'best');

                ax = nexttile(t, 2);
                hold(ax, 'on');
                local_plot_series(ax, outer_epoch, sqrt_outer_train, [0 0.447 0.741], 'train');
                local_plot_series(ax, outer_epoch, sqrt_outer_val, [0.85 0.325 0.098], 'val');
                xlabel(ax, 'outer epoch');
                ylabel(ax, 'sqrt(outer squared metric)');
                title(ax, 'Outer RMS residual');
                grid(ax, 'on');
                legend(ax, 'Location', 'best');

                ax = nexttile(t, 3);
                if ~isempty(outer_cond)
                    plot(ax, outer_epoch(1:min(numel(outer_epoch), numel(outer_cond))), ...
                        log10(max(outer_cond(1:min(numel(outer_epoch), numel(outer_cond))), realmin)), ...
                        '-o', 'LineWidth', 1.2, 'MarkerSize', 4);
                end
                xlabel(ax, 'outer epoch');
                ylabel(ax, 'log10(eigvec cond)');
                title(ax, 'Spectral conditioning');
                grid(ax, 'on');

                ax = nexttile(t, 4);
                finite_mask = isfinite(sqrt_outer_val);
                if any(finite_mask)
                    epoch_for_scatter = outer_epoch(1:numel(sqrt_outer_val));
                    scatter(ax, epoch_for_scatter(finite_mask), sqrt_outer_val(finite_mask), 24, 'filled');
                    [~, best_idx] = min(sqrt_outer_val(finite_mask));
                    finite_idx = find(finite_mask);
                    best_global_idx = finite_idx(best_idx);
                    hold(ax, 'on');
                    scatter(ax, epoch_for_scatter(best_global_idx), sqrt_outer_val(best_global_idx), 80, 'filled');
                    legend(ax, {'outer val', sprintf('best epoch %d', round(epoch_for_scatter(best_global_idx)))}, ...
                        'Location', 'best');
                end
                xlabel(ax, 'outer epoch');
                ylabel(ax, 'sqrt(val outer squared metric)');
                title(ax, 'Best validation epoch');
                grid(ax, 'on');

                diag_png = fullfile(output_root, base_tag + "__sqrt_loss_diagnostics.png");
                exportgraphics(fig, diag_png, 'Resolution', 180);
                if include_pdf
                    exportgraphics(fig, fullfile(output_root, base_tag + "__sqrt_loss_diagnostics.pdf"));
                end
                close(fig);

                fig = figure('Visible', 'off', 'Color', 'w', 'Position', [120 120 1200 760]);
                t = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
                title(t, sprintf('%s %s %s sqrt inner metrics', ...
                    dataset_stem, observable_mode, residual_form), 'Interpreter', 'none');

                ax = nexttile(t, 1);
                hold(ax, 'on');
                local_plot_series(ax, 1:numel(sqrt_inner_squared), sqrt_inner_squared, [0 0.447 0.741], 'train');
                local_plot_series(ax, 1:numel(sqrt_inner_val_squared), sqrt_inner_val_squared, [0.85 0.325 0.098], 'val');
                xlabel(ax, 'inner fit epoch');
                ylabel(ax, 'sqrt(squared)');
                title(ax, 'Inner RMS residual from squared metric');
                grid(ax, 'on');
                legend(ax, 'Location', 'best');

                ax = nexttile(t, 2);
                hold(ax, 'on');
                local_plot_series(ax, 1:numel(sqrt_inner_per_dim), sqrt_inner_per_dim, [0 0.447 0.741], 'train');
                local_plot_series(ax, 1:numel(sqrt_inner_val_per_dim), sqrt_inner_val_per_dim, [0.85 0.325 0.098], 'val');
                xlabel(ax, 'inner fit epoch');
                ylabel(ax, 'sqrt(per\_dim)');
                title(ax, 'Inner per-dimension RMS residual');
                grid(ax, 'on');
                legend(ax, 'Location', 'best');

                metrics_png = fullfile(output_root, base_tag + "__sqrt_inner_metrics.png");
                exportgraphics(fig, metrics_png, 'Resolution', 180);
                if include_pdf
                    exportgraphics(fig, fullfile(output_root, base_tag + "__sqrt_inner_metrics.pdf"));
                end
                close(fig);

                status = "exported";
                reason = "";
                export_count = export_count + 1;
            end
        end
    catch ME
        reason = string(ME.message);
    end

    if status ~= "exported"
        skip_count = skip_count + 1;
    end

    records(end + 1, :) = { ...
        status, ...
        reason, ...
        summary_path, ...
        local_guess_family(summary_path), ...
        local_guess_dataset_stem(summary_path), ...
        local_guess_observable_mode(summary_path), ...
        local_guess_residual_form(summary_path), ...
        local_take_last_path_part(string(summary_files(k).folder)), ...
        local_safe_string(loss_mode), ...
        local_safe_string(requested_loss_mode), ...
        diag_png, ...
        metrics_png ...
    };
end

record_table = cell2table(records, 'VariableNames', { ...
    'status', ...
    'reason', ...
    'summary_path', ...
    'family', ...
    'dataset_stem', ...
    'observable_mode', ...
    'residual_form', ...
    'run_label', ...
    'loss_mode', ...
    'requested_loss_mode', ...
    'sqrt_diagnostics_png', ...
    'sqrt_inner_metrics_png' ...
});
writetable(record_table, fullfile(output_root, 'sqrt_loss_figures_index.csv'));

fid = fopen(fullfile(output_root, 'README.txt'), 'w');
fprintf(fid, 'ResKoopNet sqrt loss figures (one folder)\n');
fprintf(fid, 'Created: %s\n', char(datetime('now')));
fprintf(fid, 'Source root: %s\n', source_root);
fprintf(fid, 'Exported runs: %d\n', export_count);
fprintf(fid, 'Skipped runs: %d\n', skip_count);
fprintf(fid, 'Included loss modes: squared only = %d\n', ~include_non_squared_runs);
fprintf(fid, ['Each run contributes flattened PNG/PDF exports for sqrt(loss), ' ...
    'sqrt(outer squared metric), sqrt(squared), and sqrt(per_dim).\n']);
fclose(fid);

fprintf('Saved sqrt loss figures to: %s\n', output_root);
fprintf('Exported: %d | Skipped: %d\n', export_count, skip_count);


function text_value = local_get_text_field(S, field_name, default_value)
if isstruct(S) && isfield(S, field_name)
    value = S.(field_name);
else
    value = default_value;
end
text_value = local_to_string(value, default_value);
end


function numeric_value = local_get_numeric_field(S, field_name)
numeric_value = [];
if ~isstruct(S) || ~isfield(S, field_name)
    return;
end
value = S.(field_name);
if isempty(value)
    return;
end
if iscell(value)
    if isempty(value)
        return;
    end
    value = value{1};
end
if isnumeric(value) || islogical(value)
    numeric_value = double(value(:));
end
end


function numeric_value = local_first_numeric_field(S, field_names)
numeric_value = [];
for i = 1:numel(field_names)
    numeric_value = local_get_numeric_field(S, field_names{i});
    if ~isempty(numeric_value)
        return;
    end
end
end


function out = local_sqrt_nonnegative(x)
if isempty(x)
    out = [];
    return;
end
out = x;
finite_mask = isfinite(out);
out(finite_mask) = sqrt(max(out(finite_mask), 0));
end


function local_plot_series(ax, x, y, color_value, label_text)
if isempty(x) || isempty(y)
    return;
end
n = min(numel(x), numel(y));
plot(ax, x(1:n), y(1:n), 'LineWidth', 1.2, 'Color', color_value, 'DisplayName', label_text);
end


function text_value = local_to_string(value, default_value)
if nargin < 2
    default_value = "";
end
if isempty(value)
    text_value = string(default_value);
    return;
end
if isstring(value)
    text_value = value(1);
    return;
end
if ischar(value)
    text_value = string(value);
    return;
end
if iscell(value)
    text_value = local_to_string(value{1}, default_value);
    return;
end
if isnumeric(value) || islogical(value)
    if isscalar(value)
        text_value = string(value);
    else
        text_value = string(default_value);
    end
    return;
end
text_value = string(default_value);
end


function family = local_guess_family(summary_path)
summary_path = lower(string(summary_path));
bold_token = string(filesep) + "bold" + string(filesep);
if contains(summary_path, bold_token)
    family = "BOLD";
else
    family = "BLP";
end
end


function dataset_stem = local_guess_dataset_stem(summary_path)
parts = split(string(summary_path), filesep);
dataset_stem = "unknown";
for i = 1:numel(parts)
    token = lower(parts(i));
    if any(strcmp(token, ["e10gb1", "e10fv1", "e10gh1", "f12m01"]))
        dataset_stem = parts(i);
        return;
    end
end
end


function observable_mode = local_guess_observable_mode(summary_path)
summary_path = lower(string(summary_path));
if contains(summary_path, "complex_split")
    observable_mode = "complex_split";
elseif contains(summary_path, "abs")
    observable_mode = "abs";
elseif contains(summary_path, string(filesep) + "bold" + string(filesep))
    observable_mode = "bold";
else
    observable_mode = "unknown";
end
end


function residual_form = local_guess_residual_form(summary_path)
summary_path = lower(string(summary_path));
if contains(summary_path, "projected_vlambda")
    residual_form = "projected_vlambda";
elseif contains(summary_path, "projected_kv")
    residual_form = "projected_kv";
else
    residual_form = "unknown";
end
end


function out = local_take_last_path_part(path_value)
path_text = local_to_string(path_value, "");
if strlength(path_text) == 0
    out = "";
    return;
end
parts = split(path_text, filesep);
parts = parts(parts ~= "");
if isempty(parts)
    out = path_text;
else
    out = parts(end);
end
end


function safe_name = local_sanitize_filename(value)
safe_name = local_to_string(value, "untitled");
safe_name = regexprep(safe_name, '[^\w.-]+', '_');
safe_name = regexprep(safe_name, '_+', '_');
safe_name = regexprep(safe_name, '^_+|_+$', '');
end


function out = local_safe_string(value)
out = char(local_to_string(value, ""));
end
