function S = compute_blp_region_spectrograms(D, cfg, output_root)
% Thin in-memory wrapper around the canonical streamed spectrogram pipeline.
%
% This keeps the old convenience entry point for small runs, but the real
% implementation now lives in compute_blp_region_spectrograms_streamed.
% The wrapper asks the streamed path to also load the saved arrays back into
% memory so existing callers still receive tmpall_mean_abs/tmpall_mean_complex.

if nargin < 3 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

opts = struct();
opts.return_data = true;

S = compute_blp_region_spectrograms_streamed(D, cfg, output_root, opts);

if isfield(S, 'tmpall_mean_abs') && ~isempty(S.tmpall_mean_abs)
    S.tmpall_mean_abs = double(S.tmpall_mean_abs);
end
if isfield(S, 'tmpall_mean_complex') && ~isempty(S.tmpall_mean_complex)
    S.tmpall_mean_complex = double(S.tmpall_mean_complex);
end

S.data_in_memory = true;
S.canonical_impl = 'compute_blp_region_spectrograms_streamed';
end
