% Batch correlate BOLD eigenfunctions/deconvolved eigenfunctions with BLP density.
%
% This script expects BOLD postprocessing MAT files from:
%   script_postprocess_bold_reskoopnet_results.m
%
% For each BOLD_POST file, it automatically uses the matching dataset's:
%   E:\DataPons_processed\<dataset>\event_density\<dataset>_event_density_2s.mat

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));
rehash;
clear compute_bold_efun_density_cross_correlation
clear plot_bold_efun_density_cross_correlation_summary

close all force;
set(groot, 'defaultFigureVisible', 'off');

fprintf('Using compute function:\n  %s\n', which('compute_bold_efun_density_cross_correlation'));
fprintf('Using plot function:\n  %s\n', which('plot_bold_efun_density_cross_correlation_summary'));

processed_root = get_project_processed_root();

dataset_stems = {'e10gb1', 'e10gh1', 'e10fV1', 'f12m01'};
post_folder_name = 'bold_reskoopnet_postprocessing';

% Leave empty to run every discovered BOLD_POST file.
% Example: run_name_contains = {'HP_svd100', 'projected_kv'};
run_name_contains = {};

% Recompute by default because older failed attempts may have saved empty
% placeholder MAT files.
force_recompute = true;
continue_on_error = true;
make_figures = true;

base_params = struct();
base_params.max_lag_sec = 120;
base_params.border_mask_sec = 120;
base_params.min_valid_samples = 20;
base_params.top_n = 5;
base_params.feature_names = {'efun_abs', 'efun_real', 'deconv_abs', 'deconv_real'};
base_params.save_results = true;
base_params.make_figures = make_figures;

bold_post_files = local_discover_bold_post_files( ...
    processed_root, dataset_stems, post_folder_name, run_name_contains);

if isempty(bold_post_files)
    fprintf(2, ['No BOLD postprocessing MAT files were found.\n', ...
        'Run this first:\n  run(''scripts\\script_postprocess_bold_reskoopnet_results.m'')\n']);
    error('No BOLD postprocessing MAT files found under %s.', processed_root);
end

fprintf('Discovered %d BOLD_POST file(s).\n', numel(bold_post_files));

rows = repmat(local_empty_row(), numel(bold_post_files), 1);
row_count = 0;

for i_file = 1:numel(bold_post_files)
    bold_post_file = bold_post_files{i_file};
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s\n', i_file, numel(bold_post_files), bold_post_file);

    try
        dataset_stem = local_dataset_stem_from_bold_post(bold_post_file);
        density_file = fullfile(processed_root, dataset_stem, 'event_density', ...
            sprintf('%s_event_density_2s.mat', dataset_stem));
        if exist(density_file, 'file') ~= 2
            error('Missing event density file for %s:\n  %s', dataset_stem, density_file);
        end

        save_dir = fullfile(fileparts(fileparts(bold_post_file)), ...
            'density_cross_correlation');
        save_tag = 'bold_efun_density_xcorr';
        mat_file = fullfile(save_dir, [save_tag, '.mat']);
        if ~force_recompute && exist(mat_file, 'file') == 2
            fprintf('Existing cross-correlation found; skipping.\n  %s\n', mat_file);
            row_count = row_count + 1;
            rows(row_count) = local_make_row(bold_post_file, dataset_stem, ...
                'skipped_existing', '', mat_file);
            continue;
        end

        source = struct();
        source.name = 'blp_event_density';
        source.type = 'event_density';
        source.file = density_file;
        density_sources = {source};

        params = base_params;
        params.save_dir = save_dir;
        params.save_tag = save_tag;

        out = compute_bold_efun_density_cross_correlation( ...
            bold_post_file, density_sources, params);

        row_count = row_count + 1;
        rows(row_count) = local_make_row(bold_post_file, dataset_stem, ...
            'ok', '', out.save_paths.mat_file);

        fprintf('Saved cross-correlation MAT:\n  %s\n', out.save_paths.mat_file);
        if isfield(out, 'figure_paths')
            fprintf('Summary figure:\n  %s\n', out.figure_paths.summary_png);
            fprintf('Top curves figure:\n  %s\n', out.figure_paths.top_curves_png);
            fprintf('Top overlay figure:\n  %s\n', out.figure_paths.top_overlay_png);
        end
    catch ME
        message = getReport(ME, 'extended', 'hyperlinks', 'off');
        fprintf(2, '[ERROR] %s\n', message);
        row_count = row_count + 1;
        rows(row_count) = local_make_row(bold_post_file, '', 'error', message, '');
        if ~continue_on_error
            rethrow(ME);
        end
    end
    close all force;
end

rows = rows(1:row_count);
manifest_dir = fullfile(processed_root, 'postprocessing_manifests', ...
    'bold_efun_density_cross_correlation');
if exist(manifest_dir, 'dir') ~= 7
    mkdir(manifest_dir);
end
manifest_csv = fullfile(manifest_dir, sprintf('bold_efun_density_xcorr_manifest_%s.csv', ...
    char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
writetable(struct2table(rows), manifest_csv);

fprintf('\nFinished BOLD eigenfunction-density cross-correlation batch.\n');
fprintf('Manifest CSV:\n  %s\n', manifest_csv);


function files = local_discover_bold_post_files(processed_root, dataset_stems, folder_name, run_name_contains)
files = {};
for i_ds = 1:numel(dataset_stems)
    root_dir = fullfile(processed_root, dataset_stems{i_ds}, folder_name);
    L = dir(fullfile(root_dir, '*', 'mat', '*_bold_post.mat'));
    for i = 1:numel(L)
        file = fullfile(L(i).folder, L(i).name);
        if local_matches_all(file, run_name_contains)
            files{end + 1, 1} = file; %#ok<AGROW>
        end
    end
end
end


function tf = local_matches_all(text_value, patterns)
tf = true;
for i = 1:numel(patterns)
    if ~contains(text_value, patterns{i}, 'IgnoreCase', true)
        tf = false;
        return;
    end
end
end


function dataset_stem = local_dataset_stem_from_bold_post(bold_post_file)
parts = split(string(bold_post_file), filesep);
known = ["e10gb1", "e10gh1", "e10fV1", "f12m01"];
hit = find(ismember(lower(parts), lower(known)), 1, 'last');
if ~isempty(hit)
    dataset_stem = char(parts(hit));
    return;
end

S = load(bold_post_file, 'BOLD_POST');
if isfield(S, 'BOLD_POST') && isfield(S.BOLD_POST, 'run_info') && ...
        isfield(S.BOLD_POST.run_info, 'dataset_stem') && ...
        ~isempty(S.BOLD_POST.run_info.dataset_stem)
    dataset_stem = char(string(S.BOLD_POST.run_info.dataset_stem));
    return;
end
error('Could not infer dataset_stem from %s.', bold_post_file);
end


function row = local_empty_row()
row = struct('bold_post_file', '', 'dataset_stem', '', ...
    'status', '', 'message', '', 'xcorr_mat_file', '');
end


function row = local_make_row(bold_post_file, dataset_stem, status, message, mat_file)
row = local_empty_row();
row.bold_post_file = bold_post_file;
row.dataset_stem = dataset_stem;
row.status = status;
row.message = message;
row.xcorr_mat_file = mat_file;
end
