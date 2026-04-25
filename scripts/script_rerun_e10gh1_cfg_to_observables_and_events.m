% Rerun E10.gH1 from the current cfg to observables and event outputs.
%
% Usage from MATLAB:
%   cd('D:/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished');
%   run('scripts/script_rerun_e10gh1_cfg_to_observables_and_events.m');
%
% Optional controls, set before run(...) if needed:
%   run_spectrogram = true;
%   run_dictionaries = true;
%   run_events = true;
%   force_spectrogram = true;
%   force_events = true;
%   dict_modes = {'abs', 'complex_split'};
%
% This script intentionally does not edit cfg_E10gH1.m.

clearvars -except run_spectrogram run_dictionaries run_events ...
    force_spectrogram force_events dict_modes require_expected_sessions ...
    expected_session_ids

if ~exist('run_spectrogram', 'var') || isempty(run_spectrogram)
    run_spectrogram = true;
end

if ~exist('run_dictionaries', 'var') || isempty(run_dictionaries)
    run_dictionaries = true;
end

if ~exist('run_events', 'var') || isempty(run_events)
    run_events = true;
end

if ~exist('force_spectrogram', 'var') || isempty(force_spectrogram)
    force_spectrogram = true;
end

if ~exist('force_events', 'var') || isempty(force_events)
    force_events = true;
end

if ~exist('dict_modes', 'var') || isempty(dict_modes)
    dict_modes = {'abs', 'complex_split'};
end
dict_modes = cellstr(string(dict_modes));

if ~exist('require_expected_sessions', 'var') || isempty(require_expected_sessions)
    require_expected_sessions = true;
end

if ~exist('expected_session_ids', 'var') || isempty(expected_session_ids)
    expected_session_ids = (1:14).';
end
expected_session_ids = double(expected_session_ids(:));

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

eeglab_root = 'D:\Onedrive\Toolbox\eeglab10_2_5_8b\';
if exist(eeglab_root, 'dir') == 7 && exist('finputcheck', 'file') ~= 2
    addpath(genpath(eeglab_root));
end

othercolor_root = 'D:\Onedrive\util_functions\othercolor\';
if exist(othercolor_root, 'dir') == 7
    addpath(genpath(othercolor_root));
end

cfg = cfg_E10gH1();
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';
cfg.plot.trace_scale = 0.18;
cfg.plot.trace_clip = 4;
cfg.plot.within_gap = 1.4;
cfg.plot.between_gap = 2.2;
cfg.plot.trace_linewidth = 0.4;

output_root = get_project_processed_root();
log_dir = fullfile(output_root, cfg.file_stem, 'postprocessing', 'logs');
if exist(log_dir, 'dir') ~= 7
    mkdir(log_dir);
end

run_stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
diary_file = fullfile(log_dir, sprintf('%s_rerun_cfg_observables_events_%s.log', ...
    cfg.file_stem, run_stamp));

diary(diary_file);
diary on;
diary_cleanup = onCleanup(@() diary('off'));

try
    fprintf('\n============================================================\n');
    fprintf('E10.gH1 cfg -> observables -> event rerun\n');
    fprintf('Started: %s\n', char(datetime('now')));
    fprintf('Repo root:\n  %s\n', repo_root);
    fprintf('Output root:\n  %s\n', output_root);
    fprintf('Log file:\n  %s\n', diary_file);
    fprintf('run_spectrogram=%d, run_dictionaries=%d, run_events=%d\n', ...
        run_spectrogram, run_dictionaries, run_events);
    fprintf('force_spectrogram=%d, force_events=%d\n', ...
        force_spectrogram, force_events);
    fprintf('dict_modes: %s\n', strjoin(dict_modes, ', '));

    local_require_function('timefreqMB');
    local_require_function('finputcheck');
    local_require_function('filterSignal_mirror');
    local_require_function('find_peak_loc');

    fprintf('\n--- Loading cfg-selected BLP metadata ---\n');
    load_opts = struct();
    load_opts.metadata_only = true;
    load_opts.cache_to_disk = false;
    load_opts.force_recompute = true;
    D_meta = load_blp_dataset(cfg, load_opts);
    local_print_dataset_summary(D_meta);

    if require_expected_sessions
        if ~isequal(double(D_meta.session_ids(:)), expected_session_ids)
            error(['cfg_E10gH1 currently selects sessions %s, but this rerun script ' ...
                'expects %s. Edit cfg or set require_expected_sessions=false before run.'], ...
                mat2str(double(D_meta.session_ids(:)).'), ...
                mat2str(expected_session_ids(:).'));
        end
    end

    S = [];
    if run_spectrogram
        fprintf('\n--- Stage 1: streamed region-mean spectrograms ---\n');
        spec_opts = struct();
        spec_opts.save_precision = 'single';
        spec_opts.return_data = false;
        spec_opts.force_recompute = force_spectrogram;
        S = compute_blp_region_spectrograms_streamed(D_meta, cfg, output_root, spec_opts);
        local_verify_spectrogram_files(S.abs_file, S.complex_file, D_meta.n_time);
        fprintf('Spectrogram complete.\n');
        fprintf('  ABS:     %s\n', S.abs_file);
        fprintf('  COMPLEX: %s\n', S.complex_file);
    elseif run_dictionaries
        fprintf('\n--- Stage 1 skipped: streamed region-mean spectrograms ---\n');
        fprintf('Checking existing spectrograms because dictionary building needs them.\n');
        S = local_find_expected_spectrogram_files(cfg, output_root);
        local_verify_spectrogram_files(S.abs_file, S.complex_file, D_meta.n_time);
    else
        fprintf('\n--- Stage 1 skipped: streamed region-mean spectrograms ---\n');
        if run_events
            fprintf(['Event pipeline will reuse complete spectrograms if present, ' ...
                'or recompute them if they are incomplete.\n']);
        end
    end

    dict_outputs = struct([]);
    if run_dictionaries
        fprintf('\n--- Stage 2: ResKoopNet observable dictionaries ---\n');
        dict_params = struct();
        dict_params.low_full_max_hz = 50;
        dict_params.high_max_hz = 250;
        dict_params.high_group_size = 2;
        dict_params.chunk_size = 100000;
        dict_params.precision = 'single';

        for i_mode = 1:numel(dict_modes)
            dict_params.spec_mode = dict_modes{i_mode};
            fprintf('\nBuilding dictionary mode: %s\n', dict_params.spec_mode);
            dict = build_reskoopnet_dicts(D_meta, cfg, output_root, dict_params);
            local_verify_dictionary_file(dict.save_file, D_meta.n_time, D_meta.session_ids);
            dict_outputs(i_mode).spec_mode = dict_params.spec_mode;
            dict_outputs(i_mode).save_file = dict.save_file;
            dict_outputs(i_mode).info_csv_file = dict.info_csv_file;
            dict_outputs(i_mode).n_time = dict.n_time;
            dict_outputs(i_mode).n_obs = dict.n_obs;
            fprintf('Dictionary complete: %s\n', dict.save_file);
        end
    else
        fprintf('\n--- Stage 2 skipped: ResKoopNet observable dictionaries ---\n');
    end

    clear dict

    if run_events
        fprintf('\n--- Stage 3: event detection -> consensus state diversity top windows ---\n');
        params = struct();
        params.event_params = struct( ...
            'passband', [2, 15; 30, 90; 90, 190], ...
            'band_labels', {{'theta', 'gamma', 'ripple'}}, ...
            'L_start_range', [151, 101, 51], ...
            'L_extract_range', [301, 201, 101], ...
            'ThresRatio_range', [3.5, 4, 4], ...
            'input_normalization', 'zscore_per_channel', ...
            'force_recompute', force_events);
        params.density_params = struct( ...
            'bin_sec', 2, ...
            'smooth_sigma_sec', 2, ...
            'force_recompute', force_events);
        params.consensus_params = struct( ...
            'min_channel_count', [], ...
            'require_region_presence', false, ...
            'required_regions', {{'hp', 'pl'}}, ...
            'force_recompute', force_events);
        params.summary_params = struct( ...
            'save_csv', true, ...
            'force_recompute', force_events);
        params.spectrogram_opts = struct( ...
            'save_precision', 'single', ...
            'return_data', false, ...
            'force_recompute', false);
        params.window_params = struct( ...
            'window_length_samples', 6000, ...
            'window_mode', 'global', ...
            'keep_partial_window', false, ...
            'top_k', 30, ...
            'save_csv', true, ...
            'force_recompute', force_events);
        params.plot_params = struct( ...
            'skip_existing', false);

        out = run_blp_pipeline_to_state_diversity_top_windows(cfg, output_root, params);
        local_verify_event_pipeline_output(out, D_meta.session_ids);

        fprintf('Event line complete.\n');
        fprintf('  Events:       %s\n', out.R.save_file);
        fprintf('  Density:      %s\n', out.E.save_file);
        fprintf('  Consensus:    %s\n', out.C.save_file);
        fprintf('  Summary:      %s\n', out.S.save_file);
        fprintf('  Windows:      %s\n', out.W.save_file);
        fprintf('  Top plots:    %s\n', out.P.save_dir);
        if isfield(out.W, 'top_windows_table')
            fprintf('\nTop 30 state-diversity windows:\n');
            disp(out.W.top_windows_table);
        end
    else
        fprintf('\n--- Stage 3 skipped: event detection line ---\n');
    end

    fprintf('\n============================================================\n');
    fprintf('DONE: E10.gH1 rerun completed at %s\n', char(datetime('now')));
    fprintf('Log file:\n  %s\n', diary_file);
    fprintf('============================================================\n');
catch ME
    fprintf(2, '\n============================================================\n');
    fprintf(2, 'FAILED at %s\n', char(datetime('now')));
    fprintf(2, 'Log file:\n  %s\n', diary_file);
    fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    fprintf(2, '============================================================\n');
    rethrow(ME);
end

function local_require_function(function_name)
if exist(function_name, 'file') ~= 2
    error('%s was not found on the MATLAB path.', function_name);
end
end

function local_print_dataset_summary(D)
fprintf('Selected session ids: %s\n', mat2str(double(D.session_ids(:)).'));
fprintf('Session lengths:      %s\n', mat2str(double(D.session_lengths(:)).'));
fprintf('Total n_time:         %d\n', D.n_time);
fprintf('Selected channels:    %s\n', mat2str(double(D.selected_channels(:)).'));
if isfield(D, 'session_dx') && ~isempty(D.session_dx)
    fprintf('Session dx range:     [%.12g, %.12g]\n', ...
        min(double(D.session_dx(:))), max(double(D.session_dx(:))));
end
end

function S = local_find_expected_spectrogram_files(cfg, output_root)
save_dir = fullfile(output_root, cfg.file_stem, 'spectrograms');
pad_sec = 20;
pad_mode = 'mirror';
if isfield(cfg, 'spectrogram')
    if isfield(cfg.spectrogram, 'pad_sec')
        pad_sec = cfg.spectrogram.pad_sec;
    end
    if isfield(cfg.spectrogram, 'pad_mode')
        pad_mode = cfg.spectrogram.pad_mode;
    end
end

if strcmpi(pad_mode, 'mirror')
    pad_tag = sprintf('_mirrorpad_%gs', pad_sec);
else
    pad_tag = '_nopad';
end
pad_tag = strrep(pad_tag, '.', 'p');

S = struct();
S.abs_file = fullfile(save_dir, [cfg.file_stem, pad_tag, '_regionmean_spectrograms_abs.mat']);
S.complex_file = fullfile(save_dir, [cfg.file_stem, pad_tag, '_regionmean_spectrograms_complex.mat']);
end

function local_verify_spectrogram_files(abs_file, complex_file, expected_n_time)
if exist(abs_file, 'file') ~= 2
    error('ABS spectrogram file is missing: %s', abs_file);
end
if exist(complex_file, 'file') ~= 2
    error('Complex spectrogram file is missing: %s', complex_file);
end

vars_abs = who('-file', abs_file);
vars_complex = who('-file', complex_file);
required_abs = { ...
    'tmpall_mean_abs', 'freqs', 'timesout', 'regions', ...
    'session_ids', 'session_lengths', 'session_raw_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'session_dx', 'session_fs', ...
    'pad_mode', 'pad_sec', 'tmpall_mean_abs_size'};
required_complex = { ...
    'tmpall_mean_complex', 'freqs', 'regions', ...
    'session_ids', 'session_lengths', 'session_raw_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'session_dx', 'session_fs', ...
    'pad_mode', 'pad_sec', 'tmpall_mean_complex_size'};

local_assert_required_vars(vars_abs, required_abs, abs_file);
local_assert_required_vars(vars_complex, required_complex, complex_file);

wa = whos('-file', abs_file, 'tmpall_mean_abs');
wc = whos('-file', complex_file, 'tmpall_mean_complex');
if isempty(wa) || isempty(wc)
    error('Could not inspect spectrogram arrays.');
end
if numel(wa.size) < 2 || wa.size(2) ~= expected_n_time
    error('ABS spectrogram time dimension is %d, expected %d.', wa.size(2), expected_n_time);
end
if numel(wc.size) < 2 || wc.size(2) ~= expected_n_time
    error('Complex spectrogram time dimension is %d, expected %d.', wc.size(2), expected_n_time);
end

fprintf('Verified spectrogram files. ABS size=%s, complex size=%s\n', ...
    mat2str(wa.size), mat2str(wc.size));
end

function local_assert_required_vars(actual_vars, required_vars, mat_file)
missing = setdiff(required_vars, actual_vars);
if ~isempty(missing)
    error('MAT file is incomplete: %s\nMissing variables: %s', ...
        mat_file, strjoin(missing, ', '));
end
end

function local_verify_dictionary_file(dict_file, expected_n_time, expected_session_ids)
if exist(dict_file, 'file') ~= 2
    error('Dictionary file is missing: %s', dict_file);
end

w = whos('-file', dict_file, 'obs');
if isempty(w)
    error('Dictionary file does not contain obs: %s', dict_file);
end
if w.size(1) ~= expected_n_time
    error('Dictionary obs time dimension is %d, expected %d: %s', ...
        w.size(1), expected_n_time, dict_file);
end

S = load(dict_file, 'session_ids', 'session_lengths');
if ~isfield(S, 'session_ids') || ~isequal(double(S.session_ids(:)), double(expected_session_ids(:)))
    error('Dictionary session_ids do not match the current cfg: %s', dict_file);
end
if ~isfield(S, 'session_lengths') || sum(double(S.session_lengths(:))) ~= expected_n_time
    error('Dictionary session_lengths do not sum to expected n_time: %s', dict_file);
end

fprintf('Verified dictionary file. obs size=%s\n', mat2str(w.size));
end

function local_verify_event_pipeline_output(out, expected_session_ids)
if ~isfield(out, 'R') || ~isequal(double(out.R.session_ids(:)), double(expected_session_ids(:)))
    error('Event detection output session_ids do not match the current cfg.');
end
if ~isfield(out, 'W') || ~isfield(out.W, 'top_windows_table')
    error('State-diversity output does not include top_windows_table.');
end
if height(out.W.top_windows_table) ~= 30
    error('Expected 30 top state-diversity windows, got %d.', height(out.W.top_windows_table));
end
end
