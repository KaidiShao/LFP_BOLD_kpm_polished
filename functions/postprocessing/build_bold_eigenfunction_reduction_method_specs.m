function specs = build_bold_eigenfunction_reduction_method_specs(params)
%BUILD_BOLD_EIGENFUNCTION_REDUCTION_METHOD_SPECS Canonical P9 DR methods.

if nargin < 1 || isempty(params)
    params = build_bold_eigenfunction_reduction_params();
end

specs = build_blp_eigenfunction_reduction_method_specs(params);
end
