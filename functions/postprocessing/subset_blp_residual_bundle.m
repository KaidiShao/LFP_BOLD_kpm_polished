function mode_info = subset_blp_residual_bundle(residual_bundle, max_modes)
%SUBSET_BLP_RESIDUAL_BUNDLE Reuse the preselected residual mode ordering from a workspace bundle.

bundle_info = residual_bundle.mode_info;
Ksel = min(bundle_info.n_modes, max_modes);
mode_info = bundle_info;
mode_info.selected_idx = bundle_info.selected_idx(1:Ksel);
mode_info.selected_evalues = bundle_info.selected_evalues(1:Ksel);
mode_info.lambda_d = bundle_info.lambda_d(1:Ksel);
mode_info.lambda_d_row = bundle_info.lambda_d_row(1:Ksel);
mode_info.n_modes = Ksel;
end
