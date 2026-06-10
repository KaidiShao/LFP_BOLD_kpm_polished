function existing_mat = find_existing_pipeline6_cross_correlation_result(save_dir, cfg, run_info)
%FIND_EXISTING_PIPELINE6_CROSS_CORRELATION_RESULT Find the most recent non-lagged saved result.

existing_mat = '';
run_label = local_filename_safe(run_info.run_name);
pattern = sprintf('%s_%s*_koopman_residual_edmd_top*.mat', cfg.file_stem, run_label);
L = dir(fullfile(save_dir, pattern));
if isempty(L)
    return;
end

names = {L.name};
keep = ~contains(names, '_lagged', 'IgnoreCase', true);
keep = keep & [L.bytes] > 0;
L = L(keep);
if isempty(L)
    return;
end

[~, order] = sort([L.datenum], 'descend');
L = L(order);
existing_mat = fullfile(L(1).folder, L(1).name);
end


function out = local_filename_safe(name_in)
out = regexprep(char(string(name_in)), '[^a-zA-Z0-9_\-]+', '_');
out = regexprep(out, '_+', '_');
out = regexprep(out, '^_+|_+$', '');
if isempty(out)
    out = 'untitled';
end
end
