function spec = resolve_regionmean_spectrogram_files(cfg, output_root)
% Resolve the canonical saved-file contract for region-mean spectrograms.

if nargin < 2 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

if ~isfield(cfg, 'file_stem') || isempty(cfg.file_stem)
    error('cfg.file_stem is required.');
end

pad_sec = 20;
pad_mode = 'mirror';
if isfield(cfg, 'spectrogram')
    if isfield(cfg.spectrogram, 'pad_sec') && ~isempty(cfg.spectrogram.pad_sec)
        pad_sec = cfg.spectrogram.pad_sec;
    end
    if isfield(cfg.spectrogram, 'pad_mode') && ~isempty(cfg.spectrogram.pad_mode)
        pad_mode = cfg.spectrogram.pad_mode;
    end
end

if strcmpi(pad_mode, 'mirror')
    pad_tag = sprintf('_mirrorpad_%gs', pad_sec);
else
    pad_tag = '_nopad';
end
pad_tag = strrep(pad_tag, '.', 'p');

save_dir = fullfile(output_root, cfg.file_stem, 'spectrograms');

abs_name = [cfg.file_stem, pad_tag, '_regionmean_spectrograms_abs.mat'];
complex_name = [cfg.file_stem, pad_tag, '_regionmean_spectrograms_complex.mat'];

spec = struct();
spec.pad_sec = pad_sec;
spec.pad_mode = pad_mode;
spec.pad_tag = pad_tag;
spec.save_dir = save_dir;
spec.abs_name = abs_name;
spec.complex_name = complex_name;
spec.abs_file = fullfile(save_dir, abs_name);
spec.complex_file = fullfile(save_dir, complex_name);
spec.abs_var_name = 'tmpall_mean_abs';
spec.complex_var_name = 'tmpall_mean_complex';
end
