function tf = file_signature_matches(sig_a, sig_b)
%FILE_SIGNATURE_MATCHES Compare two lightweight file signatures.

tf = isstruct(sig_a) && isstruct(sig_b);
if ~tf
    return;
end

required_fields = {'source_ref', 'exists', 'bytes', 'datenum'};
for i = 1:numel(required_fields)
    if ~isfield(sig_a, required_fields{i}) || ~isfield(sig_b, required_fields{i})
        tf = false;
        return;
    end
end

tf = strcmpi(char(string(sig_a.source_ref)), char(string(sig_b.source_ref)));
if ~tf
    return;
end

tf = logical(sig_a.exists) == logical(sig_b.exists);
if ~tf || ~logical(sig_a.exists)
    return;
end

tf = isequaln(double(sig_a.bytes), double(sig_b.bytes)) && ...
    isequaln(double(sig_a.datenum), double(sig_b.datenum));
end
