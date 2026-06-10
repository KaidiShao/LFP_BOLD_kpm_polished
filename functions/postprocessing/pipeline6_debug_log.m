function pipeline6_debug_log(params, message)
%PIPELINE6_DEBUG_LOG Append an optional debug line for pipeline 6 helpers.

if ~isfield(params, 'debug_log_file') || isempty(params.debug_log_file)
    return;
end

fid = fopen(params.debug_log_file, 'a');
if fid < 0
    return;
end
cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '[%s] %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), char(message));
end
