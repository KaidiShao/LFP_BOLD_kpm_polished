function sig = build_file_signature(source_ref)
%BUILD_FILE_SIGNATURE Build a lightweight signature for cache validation.

if nargin < 1 || isempty(source_ref)
    source_ref = '';
end

source_ref = char(string(source_ref));

sig = struct();
sig.source_ref = source_ref;
sig.exists = false;
sig.bytes = [];
sig.datenum = [];

if isempty(source_ref)
    return;
end

if exist(source_ref, 'file') == 2
    info = dir(source_ref);
    sig.exists = true;
    sig.bytes = double(info.bytes);
    sig.datenum = double(info.datenum);
end
end
