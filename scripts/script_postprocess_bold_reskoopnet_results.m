% Postprocess completed BOLD MLP ResKoopNet outputs.
%
% Outputs:
%   E:\DataPons_processed\<dataset>\pipeline7_bold_reskoopnet_postprocessing\<run>\mat
%   E:\DataPons_processed\<dataset>\pipeline7_bold_reskoopnet_postprocessing\<run>\fig
%   E:\DataPons_processed\<dataset>\pipeline7_bold_reskoopnet_postprocessing\<run>\fig\intrinsic_activation_maps

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

params = build_bold_reskoopnet_postprocessing_params();
params.autodl_roots = {'E:\autodl_results_local\bold_wsl'};
params.processed_root = io_project.get_project_processed_root();

% Leave empty to scan every completed run.
params.dataset_stems = {};
params.observable_modes = {};
params.residual_forms = {};

params.headless = true;
params.close_figures = true;

manifest = postprocess_bold_reskoopnet_results(params);

fprintf('\nFinished BOLD ResKoopNet postprocessing.\n');
fprintf('Manifest CSV:\n  %s\n', manifest.csv_file);
