function params = build_bold_eigenfunction_reduction_params()
%BUILD_BOLD_EIGENFUNCTION_REDUCTION_PARAMS Default batch params for pipeline 9.

params = struct();
params.processed_root = io_project.get_project_processed_root();

params.dataset_stems = {};
params.exclude_dataset_stems = {};
params.observable_modes = {};
params.residual_forms = {};
params.run_name_filter = {};
params.run_name_contains = {};
params.max_runs = [];

params.feature_names = {'efun_abs', 'efun_real', 'deconv_abs', 'deconv_real'};
params.method_filter = {};
params.component_count_sweep = 3:20;
params.time_component_count_sweep = [];
params.spectrum_component_count_sweep = [];

params.selection = struct();
params.selection.abs_thresh = 0.01;
params.selection.sort_by = 'preserve';
params.selection.sort_dir = 'descend';
params.selection.max_modes = Inf;

params.feature_normalization = 'maxabs_per_mode';
params.save_payload = 'compact';
params.save_v7_3 = true;
params.skip_existing = true;
params.continue_on_error = true;
params.write_manifest = true;
params.headless = true;
params.close_figures = true;

params.make_summary_plot = true;
params.plot = struct();
params.plot.max_observables = 40;
params.plot.max_source_modes = 40;
params.plot.max_components = 20;
params.plot.max_samples = 2000;
params.plot.window_start = 1;
params.plot.window_idx = [];
params.plot.select_observables_by = 'variance';
params.plot.select_modes_by = 'order';
params.plot.colormap = 'turbo';
params.plot.figure_position = [80, 80, 1280, 820];
params.plot.save_png = true;
params.plot.save_fig = false;
params.plot.resolution = 220;
end
