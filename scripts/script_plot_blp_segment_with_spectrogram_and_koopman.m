this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));
results_root = get_project_results_root(repo_root);

othercolor_root = 'D:\Onedrive\util_functions\othercolor\';
if exist(othercolor_root, 'dir') == 7
    addpath(genpath(othercolor_root));
end

close all force;

cfg = cfg_F12m01();
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';

cfg.plot.trace_scale = 0.18;
cfg.plot.trace_clip = 4;
cfg.plot.within_gap = 1.4;
cfg.plot.between_gap = 2.2;
cfg.plot.trace_linewidth = 0.4;

output_root = get_project_processed_root();
time_range_sec = [300, 320];
freq_range_to_plot = [0, 250];
color_limits = [];
event_input = 'auto';
event_colors = [];

if exist('othercolor', 'file') == 2
    spec_colormap = flipud(othercolor('Spectral10'));
else
    spec_colormap = flipud(turbo(256));
end

source_cfg = struct();
source_cfg.mode = 'chunk_dir';
source_cfg.data_dir = 'E:\autodl_results\initial_point_test1_KV';
source_cfg.edmd_file = '';
source_cfg.concat = struct();
source_cfg.concat.filename_pattern = '*_outputs_*.mat';
source_cfg.concat.variable_name = 'EDMD_outputs';
source_cfg.concat.concat_fields = {'efuns'};
source_cfg.concat.concat_dim = 1;
source_cfg.concat.required_equal_fields = {'evalues', 'kpm_modes', 'N_dict', 'residual_form'};
source_cfg.concat.allow_missing_chunks = false;
source_cfg.concat.verbose = true;
source_cfg.concat.progress_every = 50;

edmd_cfg = struct();
edmd_cfg.abs_thresh = 0.01;
edmd_cfg.sort_by = 'modulus';
edmd_cfg.sort_dir = 'descend';
edmd_cfg.max_basis = 30;
edmd_cfg.do_plot = false;

koopman_cfg = struct();
koopman_cfg.use_original_sorted = false;
koopman_cfg.mode_idx = [];
koopman_cfg.max_modes = 20;
koopman_cfg.feature = 'abs';
koopman_cfg.normalize_scope = 'window';
koopman_cfg.first_u_mode = 'phi1';
koopman_cfg.koopman_cmap = spec_colormap;
koopman_cfg.koopman_clim = [];
koopman_cfg.residual_clim = [];

save_cfg = struct();
save_cfg.do_save = false;
save_cfg.save_dir = fullfile(results_root, 'blp_segment_with_spectrogram_and_koopman');

[EDMD_outputs, ~, source_info] = load_edmd_source(source_cfg);
fprintf('Loaded EDMD source mode: %s\n', source_info.mode);
fprintf('EDMD source path:\n  %s\n', source_info.path);
fprintf('Loaded %d time samples and %d modes before sorting.\n', ...
    size(EDMD_outputs.efuns, 1), size(EDMD_outputs.efuns, 2));

[EDMD_outputs_post, ~] = postprocess_EDMD_outputs(EDMD_outputs, edmd_cfg);
koopman_cfg.EDMD_outputs = EDMD_outputs_post;

prep_cfg = struct();
prep_cfg.show_events = true;
prep_cfg.event_input = event_input;
prep_cfg.band_colors = event_colors;
prep_cfg.include_spectrogram = true;

plot_data = prepare_blp_plot_data(cfg, output_root, prep_cfg);
base_plot_cache = build_blp_plot_window_cache(plot_data, time_range_sec, freq_range_to_plot);

feature_list = {'abs', 'real'};
hfig = struct();
plot_info = struct();

for i_feature = 1:numel(feature_list)
    feature_name = feature_list{i_feature};
    feature_field = matlab.lang.makeValidName(feature_name);

    koopman_cfg_this = koopman_cfg;
    koopman_cfg_this.feature = feature_name;
    koopman_cfg_this.koopman_cmap = spec_colormap;

    if strcmpi(feature_name, 'real')
        koopman_cfg_this.koopman_clim = [-1, 1];
        koopman_cfg_this.residual_clim = [-1, 1];
    else
        koopman_cfg_this.koopman_clim = [];
        koopman_cfg_this.residual_clim = [];
    end

    [hfig.(feature_field), plot_info.(feature_field)] = plot_blp_segment_with_spectrogram_and_koopman( ...
        base_plot_cache, ...
        spec_colormap, ...
        color_limits, ...
        koopman_cfg_this);
    if ~isempty(hfig.(feature_field)) && isgraphics(hfig.(feature_field))
        set(hfig.(feature_field), ...
            'Name', sprintf('BLP segment with spectrogram and Koopman (%s)', feature_name), ...
            'NumberTitle', 'off');
    end

    if save_cfg.do_save
        if exist(save_cfg.save_dir, 'dir') ~= 7
            mkdir(save_cfg.save_dir);
        end

        file_stub = sprintf('%s_koopman_segment_%s_%s_%gs_to_%gs', ...
            cfg.file_stem, ...
            local_filename_safe(feature_name), ...
            local_filename_safe(koopman_cfg.normalize_scope), ...
            time_range_sec(1), ...
            time_range_sec(2));
        file_stub = strrep(file_stub, '.', 'p');
        exportgraphics(hfig.(feature_field), ...
            fullfile(save_cfg.save_dir, [file_stub, '.png']), 'Resolution', 300);
    end
end

disp('Plot info available in workspace variable plot_info with fields plot_info.abs and plot_info.real.');


function s = local_filename_safe(s)
s = lower(char(s));
s = strrep(s, ' ', '_');
s = strrep(s, '-', '_');
end
