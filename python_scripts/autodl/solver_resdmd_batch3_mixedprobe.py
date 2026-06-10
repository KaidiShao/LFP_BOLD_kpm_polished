"""Probe-only ResKoopNet solver variant.

This module intentionally leaves the production solver untouched.  It reuses
solver_resdmd_batch3 but disables Keras' per-inner-epoch full validation pass,
which is useful when timing large BLP/BOLD standardization probes.
"""

import numpy as np

from solver_resdmd_batch3 import KoopmanSolver as _BaseKoopmanSolver


class KoopmanSolver(_BaseKoopmanSolver):
    """Mixed-precision probe solver with no inner validation pass."""

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("training_policy", "float32")
        kwargs.setdefault("analysis_dtype", "float64")
        kwargs.setdefault("gram_dtype", "float64")
        kwargs.setdefault("spectral_dtype", "float64")
        super().__init__(*args, **kwargs)

    def train_psi(self, model, epochs, callbacks):
        """Train dictionary weights without running full validation each inner epoch."""
        callbacks = [
            cb for cb in (callbacks or [])
            if getattr(cb, "monitor", None) not in {"val_loss", "val_squared_loss"}
        ]
        history = model.fit(
            x=self.ds_train,
            epochs=epochs,
            callbacks=callbacks,
            verbose=0,
        )
        n = len(history.history.get("loss", []))
        nan_values = [np.nan] * n
        history.history.setdefault("val_loss", nan_values)
        history.history.setdefault("val_squared_loss", nan_values)
        history.history.setdefault("val_per_dim_loss", nan_values)
        history.history.setdefault("val_relative_target_loss", nan_values)
        return history
