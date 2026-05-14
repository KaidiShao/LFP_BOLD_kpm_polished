import tensorflow as tf

from solver_resdmd_batch3 import KoopmanSolver as DynamicKoopmanSolver


class KoopmanSolver(DynamicKoopmanSolver):
    """
    Engineering-oriented static-target variant of the polished TensorFlow solver.

    This keeps the current accelerated / chunked / checkpointed training
    infrastructure, but freezes the projected_vlambda training target
    (V, Lambda, and optionally K for projected_kv) after the initial
    spectral sync. In other words, inner training keeps chasing the same
    spectral target throughout the run, which is the behavior that most
    closely matches the legacy local TensorFlow solver.

    We also force the faster torch-like outer schedule (`pre_only`) so
    each outer epoch performs a single spectral refresh for evaluation /
    diagnostics, without re-syncing the training graph afterward.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.static_training_target = True
        self._static_target_synced_once = False
        # Best-checkpoint restore should recompute spectral state from the
        # restored network weights because checkpointed V/Lambda remain frozen.
        self.recompute_spectral_after_checkpoint_restore = True

    def _sync_training_spectral_state(self):
        if self.eigenvectors_var is None or self.eigenvalues_var is None:
            return

        if self._static_target_synced_once:
            return

        self.eigenvectors_var.assign(tf.cast(self.eigenvectors, self.training_complex_dtype))
        self.eigenvalues_var.assign(tf.cast(self.eigenvalues, self.training_complex_dtype))
        if self.residual_form == 'projected_kv' and self.model is not None:
            layer_k_weights = self.model.get_layer('Layer_K').weights[0]
            layer_k_weights.assign(tf.cast(self.K, layer_k_weights.dtype))
        self._static_target_synced_once = True
        print("build: frozen training spectral target at initial V/Lambda")

    def build(self, *args, **kwargs):
        self._static_target_synced_once = False
        # Keep the single-refresh outer schedule from the current torch-like
        # implementation while freezing the training target after the first sync.
        kwargs['spectral_sync_mode'] = 'pre_only'
        return super().build(*args, **kwargs)
