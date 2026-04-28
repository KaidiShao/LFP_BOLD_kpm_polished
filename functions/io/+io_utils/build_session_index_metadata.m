function [session_start_idx, session_end_idx, border_idx] = build_session_index_metadata(session_lengths)
%BUILD_SESSION_INDEX_METADATA Build concatenated-session index bookkeeping.

session_lengths = session_lengths(:);

if isempty(session_lengths)
    session_start_idx = [];
    session_end_idx = [];
    border_idx = [];
    return;
end

session_end_idx = cumsum(session_lengths);
session_start_idx = [1; session_end_idx(1:end-1) + 1];
border_idx = session_end_idx(1:end-1);
