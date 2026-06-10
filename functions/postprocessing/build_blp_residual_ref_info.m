function ref_info = build_blp_residual_ref_info(residual_bundle)
%BUILD_BLP_RESIDUAL_REF_INFO Build minimal reference metadata for workspace residual sources.

ref_info = struct();
ref_info.evalues = residual_bundle.mode_info.selected_evalues;
ref_info.kpm_modes = [];
ref_info.N_dict = [];
ref_info.residual_form = 'workspace_precomputed';
end
