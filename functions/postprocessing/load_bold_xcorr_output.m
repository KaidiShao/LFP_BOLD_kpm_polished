function out = load_bold_xcorr_output(input)
%LOAD_BOLD_XCORR_OUTPUT Normalize a saved or in-memory pipeline 8 xcorr output.

if nargin < 1 || isempty(input)
    error('xcorr input is required.');
end

if ischar(input) || isstring(input)
    file = char(string(input));
    S = load(file, 'out');
    if ~isfield(S, 'out')
        error('XCORR variable ''out'' missing in %s.', file);
    end
    out = S.out;
    out.source_file = file;
elseif isstruct(input)
    if isfield(input, 'top_table')
        out = input;
    elseif isfield(input, 'out')
        out = input.out;
    else
        error('xcorr input struct must be an xcorr result or contain field ''out''.');
    end
    if ~isfield(out, 'source_file')
        out.source_file = '';
    end
else
    error('xcorr input must be a path or struct.');
end
end
