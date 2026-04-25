this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

data_dir = 'E:\autodl_results\initial_point_test1_KV';
save_dir = data_dir;

params = struct();
params.filename_pattern = '*_outputs_*.mat';
params.variable_name = 'EDMD_outputs';
params.concat_fields = {'efuns'};
params.concat_dim = 1;
params.required_equal_fields = { ...
    'evalues', 'kpm_modes', 'N_dict', 'residual_form', ...
    'observable_tag', 'observable_mode', ...
    'dt', 'dx', 'sampling_period', 'sample_period', 'fs', 'sampling_frequency', ...
    'session_dx', 'session_fs', 'session_ids', 'session_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx'};
params.allow_missing_chunks = false;
params.verbose = true;
params.progress_every = 50;

% If future EDMD outputs also contain time-resolved fields such as Psi_X or
% Psi_Y, add them here:
% params.concat_fields = {'efuns', 'Psi_X', 'Psi_Y'};

[EDMD_outputs, concat_info] = load_and_concat_edmd_output_chunks(data_dir, params);

% if exist(save_dir, 'dir') ~= 7
%     mkdir(save_dir);
% end
% 
% save_name = sprintf('%s_outputs_%d_to_%d_concat.mat', ...
%     concat_info.file_prefix, ...
%     concat_info.chunk_ids(1), ...
%     concat_info.chunk_ids(end));
% save_path = fullfile(save_dir, save_name);
% 
% fprintf('Saving concatenated EDMD outputs to:\n  %s\n', save_path);
% save(save_path, 'EDMD_outputs', 'concat_info', '-v7.3');
% 
% fprintf('Done.\n');
% fprintf('Concatenated %d chunks into %d samples.\n', ...
%     concat_info.n_chunks, concat_info.total_length);
% fprintf('EDMD_outputs.efuns size = [%s]\n', num2str(size(EDMD_outputs.efuns)));
