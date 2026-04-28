function Spec = require_saved_blp_spectrogram_support(cfg, output_root)
%REQUIRE_SAVED_BLP_SPECTROGRAM_SUPPORT Require an existing saved region-mean
% spectrogram file for plotting support.

if nargin < 2 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

spec_files = resolve_regionmean_spectrogram_files(cfg, output_root);

abs_file = '';
complex_file = '';
if exist(spec_files.abs_file, 'file') == 2
    abs_file = spec_files.abs_file;
end
if exist(spec_files.complex_file, 'file') == 2
    complex_file = spec_files.complex_file;
end

if isempty(abs_file)
    error(['Top-window plotting requires a saved region-mean spectrogram file.\n' ...
        'Expected saved file:\n  %s\n\n' ...
        'Run the BLP spectrogram pipeline first.'], ...
        spec_files.abs_file);
end

Spec = struct();
Spec.abs_file = abs_file;
Spec.complex_file = complex_file;
Spec.pad_mode = spec_files.pad_mode;
Spec.pad_sec = spec_files.pad_sec;
Spec.source = 'saved_file_required';
end
