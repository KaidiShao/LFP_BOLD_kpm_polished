% Postprocess completed BOLD MLP ResKoopNet outputs.
%
% Outputs:
%   E:\DataPons_processed\<dataset>\bold_reskoopnet_postprocessing\<run>\mat
%   E:\DataPons_processed\<dataset>\bold_reskoopnet_postprocessing\<run>\fig

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

params = struct();
params.autodl_roots = {'E:\autodl_results_local\bold_wsl'};
params.processed_root = io_project.get_project_processed_root();

% Leave empty to scan everything that has completed output chunks.
params.dataset_stems = {};
params.observable_modes = {};
params.residual_forms = {};

params.output_folder_name = 'bold_reskoopnet_postprocessing';
params.skip_existing = true;
params.continue_on_error = true;
params.headless = true;
params.close_figures = true;

params.post_opts = struct();
params.post_opts.abs_thresh = 0.01;
params.post_opts.sort_by = 'modulus';
params.post_opts.sort_dir = 'descend';
params.post_opts.max_basis = 80;

params.deconv = struct();
params.deconv.method = 'koopman_residual';
params.deconv.lambda_source = 'edmd';
params.deconv.max_modes_all = 80;
params.deconv.max_modes_sel = 40;

params.timescale = struct();
params.timescale.max_modes_all = 80;
params.timescale.max_modes_sel = 40;
params.timescale.maxLag = 200;
params.timescale.xlim_time = 240;

manifest = postprocess_bold_reskoopnet_results(params);

fprintf('\nFinished BOLD ResKoopNet postprocessing.\n');
fprintf('Manifest CSV:\n  %s\n', manifest.csv_file);
