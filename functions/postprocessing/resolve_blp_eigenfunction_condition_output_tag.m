function condition_tag = resolve_blp_eigenfunction_condition_output_tag(run_info, params)
%RESOLVE_BLP_EIGENFUNCTION_CONDITION_OUTPUT_TAG Resolve pipeline 5 output grouping.

if nargin < 2 || isempty(params)
    params = struct();
end

if isfield(params, 'condition_output_tag') && ~isempty(params.condition_output_tag)
    condition_tag = char(string(params.condition_output_tag));
else
    mode = 'condition';
    if isfield(params, 'condition_tag_mode') && ~isempty(params.condition_tag_mode)
        mode = lower(char(string(params.condition_tag_mode)));
    end

    switch mode
        case {'condition', 'default'}
            condition_tag = build_blp_eigenfunction_condition_tag(run_info);

        case {'run', 'run_name', 'source_run'}
            condition_tag = char(string(run_info.run_name));

        case {'condition_run', 'condition_run_name'}
            condition_tag = sprintf('%s__%s', ...
                build_blp_eigenfunction_condition_tag(run_info), ...
                char(string(run_info.run_name)));

        otherwise
            error('Unsupported condition_tag_mode: %s', mode);
    end
end

condition_tag = local_filename_safe(condition_tag);
end


function out = local_filename_safe(name_in)
out = regexprep(char(string(name_in)), '[^a-zA-Z0-9_\-]+', '_');
out = regexprep(out, '_+', '_');
out = regexprep(out, '^_+|_+$', '');
if isempty(out)
    out = 'untitled_condition';
end
end
