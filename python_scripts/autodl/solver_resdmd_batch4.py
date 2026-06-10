import numpy as np
import tensorflow as tf
from tensorflow.keras import mixed_precision
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras.models import Model

import solver_resdmd_batch3 as base_solver


# Reset the global default after importing batch3, which forces float32.
tf.keras.backend.set_floatx("float64")
mixed_precision.set_global_policy("float64")


AbstractDictionary = base_solver.AbstractDictionary
KoopmanNN = base_solver.KoopmanNN


class KoopmanSolver(base_solver.KoopmanSolver):
    """
    High-precision pendulum-focused variant of solver3.

    Differences from batch3:
    - uses float64 / complex128 for training and spectral quantities
    - matches the upstream loss scaling used in solver_resdmd_tf.py
    - defaults compute_final_info() to reg_final=0.01 like upstream
    """

    def __init__(
        self,
        dic,
        target_dim,
        reg=0.01,
        training_policy="float64",
        analysis_dtype="float64",
        gram_dtype="float64",
        spectral_dtype="float64",
        residual_form="projected_vlambda",
        loss_mode="squared",
        loss_epsilon=1e-12,
    ):
        super().__init__(
            dic=dic,
            target_dim=target_dim,
            reg=reg,
            training_policy=training_policy,
            analysis_dtype=analysis_dtype,
            gram_dtype=gram_dtype,
            spectral_dtype=spectral_dtype,
            residual_form=residual_form,
            loss_mode=loss_mode,
            loss_epsilon=loss_epsilon,
        )
        self.training_input_dtype = tf.float64
        self.training_real_dtype = tf.float64
        self.training_complex_dtype = tf.complex128
        tf.keras.backend.set_floatx("float64")
        mixed_precision.set_global_policy(self.training_policy)

    def compute_final_info(self, reg_final=None):
        if reg_final is None:
            reg_final = 0.01
        return super().compute_final_info(reg_final=reg_final)

    def squared_loss(self, y_true, y_pred):
        """
        Match upstream complex_mse scaling:
        mean over batch and complex feature dimensions.
        """
        residual, _ = self._split_packed_training_output(y_pred)
        sq_sum = tf.reduce_sum(tf.square(residual), axis=-1)
        packed_dim = tf.cast(tf.shape(residual)[-1], self.training_real_dtype)
        complex_dim = tf.maximum(
            packed_dim / tf.cast(2.0, self.training_real_dtype),
            tf.cast(1.0, self.training_real_dtype),
        )
        return tf.reduce_mean(sq_sum / complex_dim)

    def per_dim_loss(self, y_true, y_pred):
        return self.squared_loss(y_true, y_pred)

    def build_model(self):
        """
        Keep solver3's packed-output plumbing, but build the graph in float64.
        """
        inputs_x = Input((self.target_dim,), dtype=self.training_input_dtype)
        inputs_y = Input((self.target_dim,), dtype=self.training_input_dtype)

        self.psi_x = self.dic_func(inputs_x)
        self.psi_y = self.dic_func(inputs_y)

        self.layer_k = None
        if self.residual_form == "projected_kv":
            self.layer_k = Dense(
                units=self.psi_y.shape[-1],
                use_bias=False,
                name="Layer_K",
                trainable=False,
                dtype=self.training_real_dtype.name,
            )

        residual, target = self._build_training_terms(self.psi_x, self.psi_y)
        packed_residual = self._pack_complex_residual(residual)
        outputs = tf.concat(
            [packed_residual, self._pack_complex_residual(target)],
            axis=-1,
        )
        return Model(inputs=[inputs_x, inputs_y], outputs=outputs)

    def _make_ds(
        self,
        x_array,
        y_array,
        shuffle=False,
        shuffle_buffer_size=None,
        shuffle_seed=None,
        reshuffle_each_iteration=True,
    ):
        x_array = np.asarray(x_array, dtype=self.training_input_dtype.as_numpy_dtype)
        y_array = np.asarray(y_array, dtype=self.training_input_dtype.as_numpy_dtype)
        num_samples = int(x_array.shape[0])
        x_dim = int(x_array.shape[1])
        y_dim = int(y_array.shape[1])

        ds = tf.data.Dataset.from_tensor_slices(np.arange(num_samples, dtype=np.int64))
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

        def map_to_zero(index_batch):
            def fetch_numpy_batches(index_values):
                return x_array[index_values], y_array[index_values]

            x_batch, y_batch = tf.numpy_function(
                fetch_numpy_batches,
                [index_batch],
                [self.training_input_dtype, self.training_input_dtype],
            )
            x_batch.set_shape([None, x_dim])
            y_batch.set_shape([None, y_dim])
            batch_size = tf.shape(index_batch)[0]
            target_dim = self.loss_output_dim or self.output_dim
            zeros = tf.zeros((batch_size, target_dim), dtype=self.training_real_dtype)
            return (x_batch, y_batch), zeros

        ds = ds.map(map_to_zero, num_parallel_calls=tf.data.AUTOTUNE)
        return ds.prefetch(tf.data.AUTOTUNE)
