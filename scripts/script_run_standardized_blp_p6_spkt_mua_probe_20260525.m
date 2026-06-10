% Run Pipeline 6 SPKT/MUA residual cross-correlation for standardized BLP runs.
%
% This branch deliberately points at the derived standardized EDMD chunk root:
%   E:\DataPons_processed\derived_autodl_results_standardize
%
% It only runs the full-time residual-channel coupling pieces requested for
% raw-vs-standardized comparison; top-window figures and timescale plots are
% left off to keep the probe focused and cheaper.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

autodl_root = 'E:\DataPons_processed\derived_autodl_results_standardize';
make_main_plot = false;
make_timescale_plot = false;
make_deconv_plot = false;
make_deconv_window_norm_plot = false;
make_spkt_residual_cross_correlation = true;
make_mua_residual_cross_correlation = true;
spkt_cross_skip_existing = false;
mua_cross_skip_existing = false;
skip_existing = true;
max_basis = 30;
continue_on_error = true; %#ok<NASGU>

run_specs = struct( ...
    'cfg_name', {'E10gb1', 'E10gH1'}, ...
    'condition_key_filter', { ...
        {{'e10gb1|abs|projected_vlambda', 'e10gb1|complex_split|projected_vlambda'}}, ...
        {{'e10gh1|abs|projected_vlambda', 'e10gh1|complex_split|projected_vlambda'}}}, ...
    'run_name_filter', { ...
        {{ ...
            'mlp_obs_blp_vlambda_abs_stdT_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allblp_e10gb1_seed1234_projected_vlambda_abs', ...
            'mlp_obs_blp_vlambda_complex_split_stdComplexPair_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allblp_e10gb1_seed1234_projected_vlambda_complex_split'}}, ...
        {{ ...
            'mlp_obs_blp_vlambda_abs_stdT_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allblp_e10gh1_seed1234_projected_vlambda_abs', ...
            'mlp_obs_blp_vlambda_complex_split_stdComplexPair_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allblp_e10gh1_seed1234_projected_vlambda_complex_split'}}} ...
    );

fprintf('Standardized BLP P6 SPKT/MUA probe\n');
fprintf('  autodl_root: %s\n', autodl_root);

for i_spec = 1:numel(run_specs)
    cfg_name = run_specs(i_spec).cfg_name; %#ok<NASGU>
    condition_key_filter = run_specs(i_spec).condition_key_filter{1}; %#ok<NASGU>
    run_name_filter = run_specs(i_spec).run_name_filter{1}; %#ok<NASGU>

    fprintf('\n[P6 standardized %d/%d] cfg=%s\n', ...
        i_spec, numel(run_specs), cfg_name);
    fprintf('  run_name_filter:\n');
    for i_run = 1:numel(run_name_filter)
        fprintf('    - %s\n', run_name_filter{i_run});
    end

    run(fullfile(repo_root, 'scripts', ...
        'script_run_one_cfg_blp_eigenfunction_postprocessing.m'));
end

fprintf('\nFinished standardized BLP P6 SPKT/MUA probe.\n');
