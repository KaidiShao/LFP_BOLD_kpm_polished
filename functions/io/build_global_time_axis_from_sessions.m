function [t, session_start_time, session_end_time, total_duration_sec] = ...
    build_global_time_axis_from_sessions(session_lengths, session_dx)
%BUILD_GLOBAL_TIME_AXIS_FROM_SESSIONS Build one continuous sample-time axis.

session_lengths = double(session_lengths(:));
session_dx = double(session_dx(:));

if numel(session_lengths) ~= numel(session_dx)
    error('session_lengths and session_dx must have the same number of elements.');
end

n_sessions = numel(session_lengths);
session_start_time = zeros(n_sessions, 1);
session_end_time = zeros(n_sessions, 1);
t_cells = cell(n_sessions, 1);
offset_sec = 0;

for k = 1:n_sessions
    n = session_lengths(k);
    dx = session_dx(k);

    if n < 0 || n ~= round(n)
        error('session_lengths must contain nonnegative integers.');
    end
    if ~isscalar(dx) || ~isfinite(dx) || dx <= 0
        error('session_dx must contain positive finite scalars.');
    end

    session_start_time(k) = offset_sec;
    if n == 0
        t_cells{k} = zeros(1, 0);
        session_end_time(k) = offset_sec;
        continue;
    end

    t_local = offset_sec + (0:n-1) * dx;
    t_cells{k} = t_local;
    session_end_time(k) = t_local(end);
    offset_sec = offset_sec + n * dx;
end

total_duration_sec = offset_sec;

if isempty(t_cells)
    t = zeros(0, 1);
else
    t = cat(2, t_cells{:});
    t = double(t(:));
end
end
