# `functions/io` Layout

This folder is organized by responsibility rather than by load order.

- `+io_raw/`
  - Raw data loaders for BLP, BOLD, and spike sessions.
- `+io_results/`
  - Loaders for saved postprocessing outputs such as events and consensus states.
- `+io_edmd/`
  - Helpers for reading and concatenating EDMD / ResKoopNet exports.
- `+io_project/`
  - Project-level path resolvers for processed/results roots.
- `+io_utils/`
  - Small I/O-adjacent helpers such as file signatures and session-axis utilities.

Top-level `.m` files in `functions/io/` are now migration stubs that error
immediately. Old unqualified calls were intentionally removed. Update call
sites to the canonical package-qualified
names in the `+io_*` subfolders, for example:

- `io_raw.load_bold_dataset(...)`
- `io_raw.load_blp_dataset(...)`
- `io_raw.load_spike_dataset(...)`
- `io_raw.read_blp_data_slice(...)`
- `io_edmd.load_edmd_source(...)`
- `io_edmd.load_and_concat_edmd_output_chunks(...)`
- `io_results.load_event_results(...)`
- `io_results.load_consensus_state_results(...)`
- `io_project.get_project_processed_root()`
- `io_project.get_project_results_root(...)`
- `io_utils.build_file_signature(...)`
- `io_utils.file_signature_matches(...)`
- `io_utils.build_global_time_axis_from_sessions(...)`
- `io_utils.resolve_global_sample_span(...)`
