function [dx, fs] = resolve_uniform_dx(session_dx, warning_message)
%RESOLVE_UNIFORM_DX Return one shared dx/fs pair when all sessions agree.

if nargin < 2 || isempty(warning_message)
    warning_message = 'Sampling period dx is inconsistent across sessions.';
end

session_dx = session_dx(:);

if isempty(session_dx)
    dx = [];
    fs = [];
    return;
end

dx0 = session_dx(1);
if all(abs(session_dx - dx0) < 1e-12)
    dx = dx0;
    fs = 1 / dx0;
else
    dx = [];
    fs = [];
    warning(warning_message);
end
