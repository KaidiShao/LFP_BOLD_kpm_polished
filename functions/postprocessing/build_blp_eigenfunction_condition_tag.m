function condition_tag = build_blp_eigenfunction_condition_tag(run_info)
%BUILD_BLP_EIGENFUNCTION_CONDITION_TAG Canonical dataset-local condition tag.

condition_tag = sprintf('%s_%s', char(run_info.observable_mode), char(run_info.residual_form));
condition_tag = regexprep(lower(condition_tag), '[^\w\-]+', '_');
end
