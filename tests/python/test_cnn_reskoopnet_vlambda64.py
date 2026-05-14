import os
import sys
import unittest
from pathlib import Path

import numpy as np


os.environ.setdefault("CUDA_VISIBLE_DEVICES", "-1")

REPO_ROOT = Path(__file__).resolve().parents[2]
CNN_ROOT = REPO_ROOT / "python_scripts" / "cnn_reskoopnet"
if str(CNN_ROOT) not in sys.path:
    sys.path.insert(0, str(CNN_ROOT))

DATA_IMPORT_ERROR = None
try:
    from cnn_reskoopnet_data_vlambda64 import (
        make_border_aware_transition_pairs,
        reshape_obs_to_spec_and_aux,
    )
except Exception as exc:  # pragma: no cover - environment dependent
    DATA_IMPORT_ERROR = exc

SOLVER_IMPORT_ERROR = None
try:
    import tensorflow as tf

    from solver_resdmd_cnn_vlambda64 import KoopmanAuxCNN, KoopmanSolver
except Exception as exc:  # pragma: no cover - environment dependent
    SOLVER_IMPORT_ERROR = exc


def make_obs_info_rows():
    return [
        {"dim_id": "1", "source": "blp", "region_label": "gb1", "part": "raw"},
        {"dim_id": "2", "source": "blp", "region_label": "gb1", "part": "raw"},
        {"dim_id": "3", "source": "spectrogram", "region_label": "gb1", "part": "hp"},
        {"dim_id": "4", "source": "spectrogram", "region_label": "gb1", "part": "hp"},
        {"dim_id": "5", "source": "spectrogram", "region_label": "gb1", "part": "pl"},
        {"dim_id": "6", "source": "spectrogram", "region_label": "gb1", "part": "pl"},
    ]


def make_obs_matrix(n_time=8):
    t = np.linspace(0.0, 1.0, n_time, dtype=np.float64)
    return np.column_stack(
        [
            np.sin(2.0 * np.pi * t),
            np.cos(2.0 * np.pi * t),
            t,
            t**2,
            np.sin(4.0 * np.pi * t),
            np.cos(4.0 * np.pi * t),
        ]
    )


@unittest.skipIf(DATA_IMPORT_ERROR is not None, f"CNN data helper dependencies unavailable: {DATA_IMPORT_ERROR}")
class TestCnnResKoopNetVlambda64(unittest.TestCase):
    def test_reshape_obs_to_spec_and_aux_respects_float64_dtype(self):
        inputs, metadata = reshape_obs_to_spec_and_aux(
            make_obs_matrix(n_time=5),
            make_obs_info_rows(),
            output_structure="dict",
            dtype=np.float64,
        )

        self.assertEqual(inputs["aux"].dtype, np.float64)
        self.assertEqual(inputs["spec"].dtype, np.float64)
        self.assertEqual(inputs["aux"].shape, (5, 2))
        self.assertEqual(inputs["spec"].shape, (5, 2, 2))
        self.assertEqual(metadata["aux_dim"], 2)
        self.assertEqual(metadata["spec_n_freq_features"], 2)
        self.assertEqual(metadata["spec_n_channels"], 2)

    @unittest.skipIf(SOLVER_IMPORT_ERROR is not None, f"TensorFlow CNN solver dependencies unavailable: {SOLVER_IMPORT_ERROR}")
    def test_solver_build_smoke_uses_float64_projected_vlambda(self):
        tf.keras.backend.clear_session()

        inputs, metadata = reshape_obs_to_spec_and_aux(
            make_obs_matrix(n_time=8),
            make_obs_info_rows(),
            output_structure="dict",
            dtype=np.float64,
        )
        data_train, data_valid, _ = make_border_aware_transition_pairs(
            inputs,
            border_idx=np.array([], dtype=int),
            train_ratio=0.6,
            seed=7,
        )

        dic = KoopmanAuxCNN(
            spec_input_shape=(metadata["spec_n_freq_features"], metadata["spec_n_channels"]),
            aux_dim=metadata["aux_dim"],
            conv_filters=(4,),
            kernel_size=1,
            pool_size=1,
            dense_units=(4,),
            n_psi_train=4,
        )
        solver = KoopmanSolver(
            dic=dic,
            target_dim=metadata["aux_dim"],
            reg=1e-3,
        )

        losses, val_losses, stop_flag, learning_rate_changes = solver.build(
            data_train=data_train,
            data_valid=data_valid,
            epochs=1,
            batch_size=2,
            lr=1e-3,
            log_interval=1,
            lr_decay_factor=0.5,
            Nepoch=1,
            end_condition=0.0,
        )

        self.assertEqual(solver.residual_form, "projected_vlambda")
        self.assertEqual(solver.training_input_dtype, tf.float64)
        self.assertEqual(solver.training_real_dtype, tf.float64)
        self.assertEqual(solver.training_complex_dtype, tf.complex128)
        self.assertGreaterEqual(len(losses), 1)
        self.assertGreaterEqual(len(val_losses), 1)
        self.assertGreaterEqual(len(learning_rate_changes), 1)
        self.assertIn(stop_flag, (0, 1))
        self.assertEqual(np.asarray(solver.eigenvalues).dtype, np.complex128)
        self.assertTrue(np.isfinite(np.asarray(losses)).all())
        self.assertTrue(np.isfinite(np.asarray(val_losses)).all())


if __name__ == "__main__":
    unittest.main()
