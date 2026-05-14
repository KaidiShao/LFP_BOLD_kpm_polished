function [EDMD_outputs, concat_info, source_info] = load_blp_eigenfunction_reduction_source(run_info, params)
%LOAD_BLP_EIGENFUNCTION_REDUCTION_SOURCE Load one completed run once for reuse across methods.

source_cfg = struct();
source_cfg.mode = 'chunk_dir';
source_cfg.data_dir = run_info.output_dir;
source_cfg.concat = struct();
source_cfg.concat.filename_pattern = params.filename_pattern;
source_cfg.concat.variable_name = params.variable_name;
source_cfg.concat.concat_fields = {'efuns'};
source_cfg.concat.concat_dim = 1;
source_cfg.concat.required_equal_fields = {'evalues', 'kpm_modes', 'N_dict', 'residual_form', ...
    'dt', 'dx', 'sampling_period', 'sample_period', 'fs', 'sampling_frequency', ...
    'session_dx', 'session_fs', 'session_ids', 'session_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx'};
source_cfg.concat.allow_missing_chunks = false;
source_cfg.concat.verbose = false;

[EDMD_outputs, concat_info, source_info] = io_edmd.load_edmd_source(source_cfg);
end
