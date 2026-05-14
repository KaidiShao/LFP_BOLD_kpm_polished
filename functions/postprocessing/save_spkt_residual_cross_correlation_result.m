function result = save_spkt_residual_cross_correlation_result(result, save_dir)
%SAVE_SPKT_RESIDUAL_CROSS_CORRELATION_RESULT Save SPKT xcorr MAT/CSV outputs.

if nargin < 1 || isempty(result)
    error('result must be provided.');
end
if nargin < 2 || isempty(save_dir)
    if isfield(result, 'params') && isfield(result.params, 'save_dir') && ~isempty(result.params.save_dir)
        save_dir = result.params.save_dir;
    else
        error('save_dir must be provided when result.params.save_dir is empty.');
    end
end

if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

source_label = local_filename_safe(local_basename(result.source.data_dir));
spike_label = local_filename_safe(result.source.spike_channel_selection_label);
tag = local_build_save_tag(result, save_dir, source_label, spike_label);

save_paths = struct();
save_paths.main_mat = fullfile(save_dir, [tag '.mat']);
save_paths.session_corr_csv = fullfile(save_dir, [tag '_session_xcorr.csv']);
save_paths.pooled_corr_csv = fullfile(save_dir, [tag '_pooled_xcorr.csv']);
save_paths.pooled_top_csv = fullfile(save_dir, [tag '_pooled_top_xcorr.csv']);
save_paths.mode_csv = fullfile(save_dir, [tag '_modes.csv']);
save_paths.session_summary_csv = fullfile(save_dir, [tag '_session_summary.csv']);

result.save_paths = save_paths;
tmp_mat = [tempname(save_dir), '.mat'];
try
    save(tmp_mat, 'result', '-v7.3');
catch ME
    if exist(tmp_mat, 'file') == 2
        delete(tmp_mat);
    end
    rethrow(ME);
end
if exist(save_paths.main_mat, 'file') == 2
    delete(save_paths.main_mat);
end
movefile(tmp_mat, save_paths.main_mat, 'f');

if ~isempty(result.session_corr_table)
    writetable(result.session_corr_table, save_paths.session_corr_csv);
else
    writetable(table(), save_paths.session_corr_csv);
end
if ~isempty(result.pooled_corr_table)
    writetable(result.pooled_corr_table, save_paths.pooled_corr_csv);
    writetable(result.pooled_top_corr, save_paths.pooled_top_csv);
else
    writetable(table(), save_paths.pooled_corr_csv);
    writetable(table(), save_paths.pooled_top_csv);
end
writetable(result.mode_table, save_paths.mode_csv);
writetable(result.session_summary, save_paths.session_summary_csv);
end


function tag = local_build_save_tag(result, save_dir, source_label, spike_label)
if contains(local_filename_safe(save_dir), source_label, 'IgnoreCase', true)
    tag = sprintf('%s_spkt_%s_top%d', ...
        result.cfg.file_stem, spike_label, height(result.mode_table));
else
    tag = sprintf('%s_%s_spkt_%s_top%d', ...
        result.cfg.file_stem, source_label, spike_label, height(result.mode_table));
end
tag = local_filename_safe(tag);
if strlength(fullfile(save_dir, [tag '.mat'])) > 240
    tag = local_filename_safe(sprintf('%s_spkt_%s_top%d', ...
        result.cfg.file_stem, spike_label, height(result.mode_table)));
end
end


function out = local_filename_safe(name_in)
out = regexprep(char(name_in), '[^a-zA-Z0-9_\-]+', '_');
out = regexprep(out, '_+', '_');
out = regexprep(out, '^_+|_+$', '');
if isempty(out)
    out = 'untitled';
end
end


function name = local_basename(path_in)
[~, name, ext] = fileparts(path_in);
name = [name, ext];
end
