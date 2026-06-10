function method_tag = build_blp_eigenfunction_reduction_method_tag(spec)
%BUILD_BLP_EIGENFUNCTION_REDUCTION_METHOD_TAG Canonical short tag for one DR method.

switch lower(char(spec.method))
    case {'svd', 'logsvd', 'nmf', 'mds', 'umap'}
        method_tag = lower(char(spec.method));
    otherwise
        method_tag = sprintf('%s_%s', lower(char(spec.kind)), lower(char(spec.method)));
end

if isfield(spec, 'options') && isstruct(spec.options) && ...
        isfield(spec.options, 'n_components') && ~isempty(spec.options.n_components)
    k = double(spec.options.n_components);
    if isfinite(k) && k > 0
        method_tag = sprintf('%s_k%02d', method_tag, round(k));
    end
end

method_tag = regexprep(method_tag, '[^\w\-]+', '_');
end
