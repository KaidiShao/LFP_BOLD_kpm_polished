function result = load_saved_pipeline6_cross_correlation_result(main_mat)
%LOAD_SAVED_PIPELINE6_CROSS_CORRELATION_RESULT Load saved result and restore standard save_paths fields.

S = load(main_mat, 'result');
if ~isfield(S, 'result')
    error('Saved cross-correlation MAT lacks variable result: %s', main_mat);
end

result = S.result;
if isfield(result, 'save_paths') && ~isempty(result.save_paths)
    return;
end

[save_dir, stem] = fileparts(main_mat);
result.save_paths = struct();
result.save_paths.main_mat = main_mat;
result.save_paths.session_corr_csv = fullfile(save_dir, [stem, '_session_xcorr.csv']);
result.save_paths.pooled_corr_csv = fullfile(save_dir, [stem, '_pooled_xcorr.csv']);
result.save_paths.pooled_top_csv = fullfile(save_dir, [stem, '_pooled_top_xcorr.csv']);
result.save_paths.mode_csv = fullfile(save_dir, [stem, '_modes.csv']);
result.save_paths.session_summary_csv = fullfile(save_dir, [stem, '_session_summary.csv']);
end
