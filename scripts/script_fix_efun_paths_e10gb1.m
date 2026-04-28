% Fix saved e10gb1 eigenfunction output paths after shortening folders/files.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

dataset_name = 'e10gb1';
output_run_name = 'mlp260415_kvabs';
new_root = fullfile(io_project.get_project_processed_root(), dataset_name, ...
    'efun', output_run_name);

result_file = local_latest_mat(fullfile(new_root, 'mat'), ...
    sprintf('%s_efun__*.mat', dataset_name));
peak_file = fullfile(new_root, 'peaks', sprintf('%s_peaks.mat', dataset_name));

fprintf('Updating result MAT paths:\n  %s\n', result_file);
S = load(result_file, 'result');
result = local_shorten_value(S.result, new_root);
result = local_set_result_paths(result, new_root, result_file, output_run_name);
save(result_file, 'result', '-v7.3');

if exist(peak_file, 'file') == 2
    fprintf('Updating peak-analysis MAT paths:\n  %s\n', peak_file);
    S = load(peak_file, 'A');
    A = local_shorten_value(S.A, new_root);
    A.source.result_file = result_file;
    A.params.save_dir = fullfile(new_root, 'peaks');
    A.params.save_tag = sprintf('%s_peaks', dataset_name);
    A.save_paths = local_peak_save_paths(new_root, dataset_name);
    save(peak_file, 'A', '-v7.3');
end

fprintf('Short-path migration complete:\n  %s\n', new_root);


function result = local_set_result_paths(result, new_root, result_file, output_run_name)
result.cfg.output.root = new_root;
result.cfg.output.figure_dir = fullfile(new_root, 'fig');
result.cfg.output.result_dir = fullfile(new_root, 'mat');
result.cfg.output.output_run_name = output_run_name;

result.cfg.save.dir = fullfile(new_root, 'mat');
result.cfg.save.file_stem = 'e10gb1_efun';
result.cfg.save.tag = 'min';

result.artifacts.result_mat_file = result_file;

if isfield(result.cfg, 'viz')
    result.cfg.viz = local_set_viz_save(result.cfg.viz, ...
        'overview', fullfile(new_root, 'fig'), 'ov');
    result.cfg.viz = local_set_viz_save(result.cfg.viz, ...
        'spectrum', fullfile(new_root, 'fig'), 'spec');
    result.cfg.viz = local_set_viz_save(result.cfg.viz, ...
        'state_space', fullfile(new_root, 'fig'), 'ss');
    result.cfg.viz = local_set_viz_save(result.cfg.viz, ...
        'state_space_consensus', fullfile(new_root, 'fig'), 'ssc');
end
end


function viz = local_set_viz_save(viz, name, save_dir, save_tag)
if isfield(viz, name) && isstruct(viz.(name))
    viz.(name).save_dir = save_dir;
    viz.(name).save_tag = save_tag;
end
end


function save_paths = local_peak_save_paths(new_root, dataset_name)
tag = sprintf('%s_peaks', dataset_name);
save_dir = fullfile(new_root, 'peaks');
save_paths = struct();
save_paths.main_mat = fullfile(save_dir, [tag, '.mat']);
save_paths.window_csv = fullfile(save_dir, [tag, '_win.csv']);
save_paths.event_peak_csv = fullfile(save_dir, [tag, '_event.csv']);
save_paths.baseline_peak_csv = fullfile(save_dir, [tag, '_base.csv']);
save_paths.stats_csv = fullfile(save_dir, [tag, '_stats.csv']);
save_paths.peak_distribution_png = fullfile(save_dir, [tag, '_dist.png']);
save_paths.mean_peak_heatmap_png = fullfile(save_dir, [tag, '_mean.png']);
save_paths.baseline_effect_heatmap_png = fullfile(save_dir, [tag, '_effect.png']);
end


function mat_file = local_latest_mat(mat_dir, pattern)
if exist(mat_dir, 'dir') ~= 7
    error('MAT directory does not exist:\n  %s', mat_dir);
end

L = dir(fullfile(mat_dir, pattern));
if isempty(L)
    error('No MAT files matching %s were found in:\n  %s', pattern, mat_dir);
end

[~, idx] = max([L.datenum]);
mat_file = fullfile(L(idx).folder, L(idx).name);
end


function value = local_shorten_value(value, new_root)
if isstruct(value)
    fields = fieldnames(value);
    for i = 1:numel(value)
        for j = 1:numel(fields)
            f = fields{j};
            value(i).(f) = local_shorten_value(value(i).(f), new_root);
        end
    end
elseif istable(value)
    names = value.Properties.VariableNames;
    for i = 1:numel(names)
        value.(names{i}) = local_shorten_value(value.(names{i}), new_root);
    end
elseif iscell(value)
    for i = 1:numel(value)
        value{i} = local_shorten_value(value{i}, new_root);
    end
elseif isstring(value)
    for i = 1:numel(value)
        value(i) = string(local_shorten_text(char(value(i)), new_root));
    end
elseif ischar(value)
    value = local_shorten_text(value, new_root);
end
end


function text = local_shorten_text(text, new_root)
old_roots = { ...
    ['E:\DataPons_processed\e10gb1\eigenfunction_reduction\', ...
    'mlp_obs_e10gb1_260415_shuffle_plateau_projected_kv_abs'], ...
    ['E:\DataPons_processed\e10gb1\koopman_postprocessing\', ...
    'eigenfunction_reduction\', ...
    'mlp_obs_e10gb1_260415_shuffle_plateau_projected_kv_abs']};

for i = 1:numel(old_roots)
    text = strrep(text, old_roots{i}, new_root);
end

text = strrep(text, '\figures\', '\fig\');
text = strrep(text, '\component_peaks_by_consensus_state\', '\peaks\');
text = strrep(text, '\top_consensus_state_diversity_eigenfunction_windows\', '\top30\');

text = regexprep(text, 'rank_(\d{2})_globalwin_(\d{3})_state_space_consensus', ...
    'r$1_w$2_ssc');
text = regexprep(text, 'rank_(\d{2})_globalwin_(\d{3})_state_space', ...
    'r$1_w$2_ss');
text = regexprep(text, 'rank_(\d{2})_globalwin_(\d{3})_overview', ...
    'r$1_w$2_ov');

text = strrep(text, '__state_space_consensus__', '__ssc__');
text = strrep(text, '__state_space__', '__ss__');
text = strrep(text, '__overview__', '__ov__');
text = strrep(text, '__spectrum_diag__', '__spec__');

text = strrep(text, '_peak_distributions', '_dist');
text = strrep(text, '_mean_peak_heatmap', '_mean');
text = strrep(text, '_baseline_effect_heatmap', '_effect');
text = strrep(text, '_event_peaks_long', '_event');
text = strrep(text, '_baseline_peaks_long', '_base');
end
