function result = build_blp_raw_eigenfunction_threshold_result(prep, run_info)
%BUILD_BLP_RAW_EIGENFUNCTION_THRESHOLD_RESULT Build minimal run-level result stub.

condition_tag = build_blp_eigenfunction_condition_tag(run_info);

result = struct();
result.meta = struct();
result.meta.method = condition_tag;

result.input = struct();
result.input.mode_index = (1:numel(prep.evalues_discrete)).';
end
