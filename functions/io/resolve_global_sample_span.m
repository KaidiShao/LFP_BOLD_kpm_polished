function span = resolve_global_sample_span( ...
    global_start, global_end, session_ids, session_lengths, session_dx, session_start_idx)
%RESOLVE_GLOBAL_SAMPLE_SPAN Map a global sample span onto session metadata.

global_start = double(global_start);
global_end = double(global_end);
session_ids = double(session_ids(:));
session_lengths = double(session_lengths(:));
session_dx = double(session_dx(:));
session_start_idx = double(session_start_idx(:));
session_end_idx = session_start_idx + session_lengths - 1;

if global_end < global_start
    error('global_end must be greater than or equal to global_start.');
end

start_k = find(global_start >= session_start_idx & global_start <= session_end_idx, 1, 'first');
end_k = find(global_end >= session_start_idx & global_end <= session_end_idx, 1, 'first');

if isempty(start_k) || isempty(end_k)
    error('Could not resolve session span for global window [%d, %d].', global_start, global_end);
end

duration_sec = 0;
for k = start_k:end_k
    overlap_start = max(global_start, session_start_idx(k));
    overlap_end = min(global_end, session_end_idx(k));
    if overlap_end >= overlap_start
        duration_sec = duration_sec + (overlap_end - overlap_start + 1) * session_dx(k);
    end
end

span = struct();
span.start_session_idx = start_k;
span.end_session_idx = end_k;
span.start_session_id = session_ids(start_k);
span.end_session_id = session_ids(end_k);
span.start_local_idx = global_start - session_start_idx(start_k) + 1;
span.end_local_idx = global_end - session_start_idx(end_k) + 1;
span.crosses_session_boundary = start_k ~= end_k;
span.duration_sec = duration_sec;
end
