"""GPU-resident data probe solver.

This solver is intentionally separate from the production solver.  It keeps the
same optimization surface as solver_resdmd_batch3_mixedprobe, but when snapshot
pairs are IndexedRowsView objects it caches one float32 copy of the base data as
a TensorFlow tensor and gathers x/y batches by row index.  That avoids repeatedly
materializing NumPy pair batches and copying them into TensorFlow.
"""

import time

import numpy as np
import tensorflow as tf
from tqdm import tqdm

from solver_resdmd_batch3_mixedprobe import KoopmanSolver as _BaseKoopmanSolver


class KoopmanSolver(_BaseKoopmanSolver):
    """Probe solver for GPU-resident IndexedRowsView data."""

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("training_policy", "float32")
        kwargs.setdefault("analysis_dtype", "float64")
        kwargs.setdefault("gram_dtype", "float64")
        kwargs.setdefault("spectral_dtype", "float64")
        super().__init__(*args, **kwargs)
        self._gpu_base_tensors = {}
        self._gpu_row_tensors = {}

    @staticmethod
    def _is_indexed_view(value):
        return (
            hasattr(value, "get_rows")
            and hasattr(value, "row_indices")
            and hasattr(value, "data")
            and hasattr(value, "shape")
        )

    def _base_tensor_for_view(self, view):
        key = id(view.data)
        tensor = self._gpu_base_tensors.get(key)
        if tensor is None:
            arr = np.asarray(view.data, dtype=self.training_input_dtype.as_numpy_dtype)
            with tf.device("/GPU:0"):
                tensor = tf.constant(arr, dtype=self.training_input_dtype)
            self._gpu_base_tensors[key] = tensor
            print(
                "mixedgpu: cached base data tensor "
                f"shape={tensor.shape}, dtype={tensor.dtype.name}, "
                f"size={arr.nbytes / (1024 ** 3):.2f} GiB"
            )
        return tensor

    def _row_tensor_for_view(self, view):
        key = id(view)
        tensor = self._gpu_row_tensors.get(key)
        if tensor is None:
            with tf.device("/GPU:0"):
                tensor = tf.constant(np.asarray(view.row_indices, dtype=np.int64), dtype=tf.int64)
            self._gpu_row_tensors[key] = tensor
        return tensor

    def _can_use_gpu_view(self, x_array, y_array):
        return (
            self._is_indexed_view(x_array)
            and self._is_indexed_view(y_array)
            and x_array.data is y_array.data
        )

    def _make_ds(
        self,
        x_array,
        y_array,
        shuffle=False,
        shuffle_buffer_size=None,
        shuffle_seed=None,
        reshuffle_each_iteration=True,
    ):
        if not self._can_use_gpu_view(x_array, y_array):
            return super()._make_ds(
                x_array,
                y_array,
                shuffle=shuffle,
                shuffle_buffer_size=shuffle_buffer_size,
                shuffle_seed=shuffle_seed,
                reshuffle_each_iteration=reshuffle_each_iteration,
            )

        base = self._base_tensor_for_view(x_array)
        x_rows = self._row_tensor_for_view(x_array)
        y_rows = self._row_tensor_for_view(y_array)
        num_samples = int(x_array.shape[0])
        x_dim = int(x_array.shape[1])
        y_dim = int(y_array.shape[1])

        ds = tf.data.Dataset.from_tensor_slices((x_rows, y_rows))
        if shuffle:
            effective_buffer_size = shuffle_buffer_size
            if effective_buffer_size is None:
                effective_buffer_size = num_samples
            effective_buffer_size = max(1, min(int(effective_buffer_size), num_samples))
            ds = ds.shuffle(
                effective_buffer_size,
                seed=shuffle_seed,
                reshuffle_each_iteration=reshuffle_each_iteration,
            )
        ds = ds.batch(self.batch_size)

        def map_to_zero(x_idx, y_idx):
            x_batch = tf.gather(base, x_idx)
            y_batch = tf.gather(base, y_idx)
            x_batch.set_shape([None, x_dim])
            y_batch.set_shape([None, y_dim])
            target_dim = self.loss_output_dim or self.output_dim
            zeros = tf.zeros((tf.shape(x_idx)[0], target_dim), dtype=self.training_real_dtype)
            return (x_batch, y_batch), zeros

        ds = ds.map(map_to_zero, num_parallel_calls=tf.data.AUTOTUNE)
        return ds.prefetch(tf.data.AUTOTUNE)

    def _view_batch(self, view, start, end, dtype):
        if not self._is_indexed_view(view):
            return tf.convert_to_tensor(view[start:end], dtype=dtype)
        base = self._base_tensor_for_view(view)
        rows = self._row_tensor_for_view(view)
        batch_idx = rows[start:end]
        return tf.cast(tf.gather(base, batch_idx), dtype)

    def compute_K(self, dic, data_x, data_y, reg):
        if not self._can_use_gpu_view(data_x, data_y):
            return super().compute_K(dic, data_x, data_y, reg)

        compute_k_start = time.time()
        sample = self._view_batch(data_x, 0, 1, self.training_input_dtype)
        analysis_dic = self._refresh_analysis_dictionary(sample)
        d = analysis_dic(tf.cast(sample, self.analysis_dtype)).shape[-1]

        XtX_sum = tf.zeros((d, d), dtype=self.gram_dtype)
        XtY_sum = tf.zeros((d, d), dtype=self.gram_dtype)
        n_total = 0
        total = int(data_x.shape[0])
        print(self.batch_size, total)

        for start in tqdm(range(0, total, self.batch_size), desc="compute_K batches", unit="batch"):
            end = min(start + self.batch_size, total)
            Xb = self._view_batch(data_x, start, end, self.analysis_dtype)
            Yb = self._view_batch(data_y, start, end, self.analysis_dtype)
            psi_x_b = tf.cast(analysis_dic(Xb), self.gram_dtype)
            psi_y_b = tf.cast(analysis_dic(Yb), self.gram_dtype)
            XtX_sum += tf.matmul(psi_x_b, psi_x_b, transpose_a=True)
            XtY_sum += tf.matmul(psi_x_b, psi_y_b, transpose_a=True)
            n_total += int(end - start)

        if n_total == 0:
            raise ValueError("compute_K received an empty dataset.")

        G = XtX_sum
        A = XtY_sum
        xtx_diag = tf.linalg.diag_part(G)
        xtx_diag_min = float(tf.reduce_min(xtx_diag).numpy())
        xtx_diag_max = float(tf.reduce_max(xtx_diag).numpy())
        xtx_has_nan = bool(tf.reduce_any(tf.math.is_nan(G)).numpy())
        xtx_has_inf = bool(tf.reduce_any(tf.math.is_inf(G)).numpy())
        a_has_nan = bool(tf.reduce_any(tf.math.is_nan(A)).numpy())
        a_has_inf = bool(tf.reduce_any(tf.math.is_inf(A)).numpy())
        print(
            "compute_K: batches done, "
            f"d={d}, n_total={n_total}, G diag min={xtx_diag_min:.3e}, G diag max={xtx_diag_max:.3e}, "
            f"G has_nan={xtx_has_nan}, G has_inf={xtx_has_inf}, "
            f"A has_nan={a_has_nan}, A has_inf={a_has_inf}, "
            f"elapsed={time.time() - compute_k_start:.2f}s"
        )

        idmat = tf.eye(d, dtype=self.gram_dtype)
        system = G + tf.cast(reg, self.gram_dtype) * idmat
        system_has_nan = bool(tf.reduce_any(tf.math.is_nan(system)).numpy())
        system_has_inf = bool(tf.reduce_any(tf.math.is_inf(system)).numpy())
        if xtx_has_nan or xtx_has_inf or a_has_nan or a_has_inf or system_has_nan or system_has_inf:
            raise FloatingPointError(
                "compute_K produced non-finite normal-equation terms. "
                f"G(has_nan={xtx_has_nan}, has_inf={xtx_has_inf}), "
                f"A(has_nan={a_has_nan}, has_inf={a_has_inf}), "
                f"system(has_nan={system_has_nan}, has_inf={system_has_inf})"
            )

        solve_start = time.time()
        print("compute_K: starting pseudoinverse solve")
        inv_term = tf.linalg.pinv(system)
        K_reg = tf.matmul(inv_term, A)
        print(
            "compute_K: pseudoinverse solve done "
            f"in {time.time() - solve_start:.2f}s, "
            f"total compute_K time={time.time() - compute_k_start:.2f}s"
        )
        k_has_nan = bool(tf.reduce_any(tf.math.is_nan(K_reg)).numpy())
        k_has_inf = bool(tf.reduce_any(tf.math.is_inf(K_reg)).numpy())
        if k_has_nan or k_has_inf:
            raise FloatingPointError(
                "compute_K produced a non-finite Koopman matrix. "
                f"has_nan={k_has_nan}, has_inf={k_has_inf}"
            )
        return K_reg

    def evaluate_spectral_losses(self, data_x, data_y, residual_form=None):
        if not self._can_use_gpu_view(data_x, data_y):
            return super().evaluate_spectral_losses(data_x, data_y, residual_form=residual_form)

        residual_form = residual_form or self.residual_form
        total_sq = 0.0
        total_target_sq = 0.0
        total_n = 0
        residual_dim = None
        batch_size = getattr(self, "batch_size", 1024)
        eigvecs = tf.cast(self.eigenvectors, self.training_complex_dtype)
        eigvals = tf.cast(self.eigenvalues, self.training_complex_dtype)
        k_matrix = tf.cast(self.K, self.training_real_dtype)

        for start in range(0, int(data_x.shape[0]), batch_size):
            end = min(start + batch_size, int(data_x.shape[0]))
            x_batch = self._view_batch(data_x, start, end, self.training_input_dtype)
            y_batch = self._view_batch(data_y, start, end, self.training_input_dtype)
            psi_x = self.dic_func(x_batch)
            psi_y = self.dic_func(y_batch)
            psi_y_v = tf.matmul(tf.cast(psi_y, self.training_complex_dtype), eigvecs)

            if residual_form == "projected_kv":
                psi_next = tf.matmul(tf.cast(psi_x, self.training_real_dtype), k_matrix)
                residual = tf.matmul(
                    tf.cast(tf.cast(psi_y, self.training_real_dtype) - psi_next, self.training_complex_dtype),
                    eigvecs,
                )
            else:
                psi_x_v = tf.matmul(tf.cast(psi_x, self.training_complex_dtype), eigvecs)
                residual = psi_y_v - psi_x_v * eigvals

            packed = self._pack_complex_residual(residual)
            target_packed = self._pack_complex_residual(psi_y_v)
            if residual_dim is None:
                residual_dim = max(1, int(packed.shape[-1]) // 2)
            batch_sq = tf.reduce_sum(tf.square(packed), axis=-1)
            target_sq = tf.reduce_sum(tf.square(target_packed), axis=-1)
            total_sq += float(tf.reduce_sum(batch_sq).numpy())
            total_target_sq += float(tf.reduce_sum(target_sq).numpy())
            total_n += int(end - start)

        if total_n == 0:
            raise ValueError("Cannot evaluate losses on an empty dataset.")

        squared = total_sq / float(total_n)
        per_dim = squared / float(residual_dim or 1)
        relative = squared / max(total_target_sq / float(total_n), self.loss_epsilon)
        return {
            "squared": float(squared),
            "per_dim": float(per_dim),
            "relative_target": float(relative),
        }
