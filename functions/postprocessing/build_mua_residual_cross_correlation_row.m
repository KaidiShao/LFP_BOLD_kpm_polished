function mua_row = build_mua_residual_cross_correlation_row(run_info, save_dir, result, overview_png, overview_fig)
%BUILD_MUA_RESIDUAL_CROSS_CORRELATION_ROW Build one MUA result summary row.

out = struct();
out.dataset = string(run_info.dataset);
out.observable_mode = string(run_info.observable_mode);
out.residual_form = string(run_info.residual_form);
out.run_name = string(run_info.run_name);
out.run_output_dir = string(run_info.output_dir);
out.mua_cross_dir = string(save_dir);
out.mua_result_mat = string(result.save_paths.main_mat);
out.pooled_xcorr_csv = string(result.save_paths.pooled_corr_csv);
out.pooled_top_xcorr_csv = string(result.save_paths.pooled_top_csv);
out.session_xcorr_csv = string(result.save_paths.session_corr_csv);
out.mode_csv = string(result.save_paths.mode_csv);
out.session_summary_csv = string(result.save_paths.session_summary_csv);
out.overview_png = string(overview_png);
out.overview_fig = string(overview_fig);
out.n_pooled_xcorr_rows = height(result.pooled_corr_table);
out.n_pooled_top_xcorr_rows = height(result.pooled_top_corr);
out.top_abs_corr = NaN;
out.top_corr = NaN;
out.top_channel_label = "";
out.top_pairing_label = "";
out.top_mode_rank = NaN;
out.status = "ok";
out.error_message = "";
if ~isempty(result.pooled_top_corr)
    top_row = result.pooled_top_corr(1, :);
    out.top_abs_corr = double(top_row.abs_corr);
    out.top_corr = double(top_row.corr);
    out.top_channel_label = string(top_row.channel_label);
    out.top_pairing_label = string(top_row.pairing_label);
    out.top_mode_rank = double(top_row.mode_rank);
end
mua_row = struct2table(out, 'AsArray', true);
end
