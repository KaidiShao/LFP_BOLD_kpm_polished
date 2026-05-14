function [EDMD_full, concat_info] = load_blp_run_edmd_outputs(run_info, params)
%LOAD_BLP_RUN_EDMD_OUTPUTS Load one completed run into a single EDMD struct.

concat_params = struct();
concat_params.filename_pattern = params.filename_pattern;
concat_params.variable_name = params.variable_name;
concat_params.concat_fields = {'efuns'};
concat_params.verbose = false;

[EDMD_full, concat_info] = io_edmd.load_and_concat_edmd_output_chunks( ...
    run_info.output_dir, concat_params);
end
