function specs = build_blp_eigenfunction_reduction_method_specs(params)
%BUILD_BLP_EIGENFUNCTION_REDUCTION_METHOD_SPECS Canonical method list for pipeline 5.

if nargin < 1
    params = struct();
end

base_specs = repmat(local_empty_method_spec(), 5, 1);
base_specs(1) = local_make_method_spec(true, 'time', 'SVD', struct());
base_specs(2) = local_make_method_spec(true, 'time', 'logSVD', struct());
base_specs(3) = local_make_method_spec(true, 'time', 'NMF', struct());
base_specs(4) = local_make_method_spec(true, 'spectrum', 'MDS', struct());
base_specs(5) = local_make_method_spec(true, 'spectrum', 'UMAP', struct());

specs = local_expand_component_counts(base_specs, params);

if isfield(params, 'method_filter') && ~isempty(params.method_filter)
    keep = false(numel(specs), 1);
    for i = 1:numel(specs)
        keep(i) = local_matches_method_filter(specs(i), params.method_filter);
    end
    specs = specs(keep);
end
end


function spec = local_empty_method_spec()
spec = struct();
spec.enabled = false;
spec.kind = '';
spec.method = '';
spec.options = struct();
end


function spec = local_make_method_spec(enabled, kind, method, options)
spec = local_empty_method_spec();
spec.enabled = enabled;
spec.kind = char(kind);
spec.method = char(method);
spec.options = options;
end


function specs = local_expand_component_counts(base_specs, params)
specs = repmat(local_empty_method_spec(), 0, 1);

for i = 1:numel(base_specs)
    base_spec = base_specs(i);
    counts = local_component_counts_for_kind(params, base_spec.kind);
    for j = 1:numel(counts)
        spec = base_spec;
        spec.options.n_components = counts(j);
        specs(end + 1, 1) = spec; %#ok<AGROW>
    end
end
end


function counts = local_component_counts_for_kind(params, kind)
counts = [];

switch lower(char(kind))
    case 'time'
        if isfield(params, 'time_component_count_sweep') && ...
                ~isempty(params.time_component_count_sweep)
            counts = params.time_component_count_sweep;
        end

    case 'spectrum'
        if isfield(params, 'spectrum_component_count_sweep') && ...
                ~isempty(params.spectrum_component_count_sweep)
            counts = params.spectrum_component_count_sweep;
        end
end

if isempty(counts) && isfield(params, 'component_count_sweep') && ...
        ~isempty(params.component_count_sweep)
    counts = params.component_count_sweep;
end

counts = local_normalize_component_counts(counts);
if isempty(counts)
    counts = 4;
end
end


function counts = local_normalize_component_counts(counts_in)
counts = double(counts_in(:).');
counts = counts(isfinite(counts) & counts > 0);
counts = unique(round(counts), 'stable');
end


function tf = local_matches_method_filter(spec, filters_in)
filters = cellstr(string(filters_in(:)));
expanded_tag = build_blp_eigenfunction_reduction_method_tag(spec);
base_tag = local_base_method_tag(spec);
candidates = unique([{expanded_tag}; {base_tag}; {lower(char(spec.method))}], 'stable');

tf = false;
for i = 1:numel(filters)
    filter_i = lower(char(string(filters{i})));
    if any(strcmpi(filter_i, candidates))
        tf = true;
        return;
    end
end
end


function tag = local_base_method_tag(spec)
tag = lower(char(spec.method));
tag = regexprep(tag, '[^\w\-]+', '_');
end
