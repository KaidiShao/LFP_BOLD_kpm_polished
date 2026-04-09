import os
# from autograd import jacobian, hessian
import time
import tensorflow as tf
from tensorflow.keras.layers import (Dense, Layer, Concatenate, Input, Conv1D, Conv2D, Conv3D,
                                     MaxPool1D, MaxPool2D, MaxPool3D,
                                     GlobalAveragePooling1D, GlobalAveragePooling2D,
                                     GlobalAveragePooling3D, BatchNormalization,
                                     Activation, Flatten)
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import mixed_precision
from tensorflow.python.ops.numpy_ops import np_config
from memory_profiler import memory_usage
import gc
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.models import load_model
import numpy as np
import scipy.linalg
from tqdm import tqdm


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.keras.backend.set_floatx('float32')

DEFAULT_AUTODL_TMP_ROOT = '/root/autodl-tmp'
DEFAULT_CHECKPOINT_ROOT = os.path.join(DEFAULT_AUTODL_TMP_ROOT, 'checkpoints')
READ_ONLY_AUTODL_ROOT = '/root/autodl-pub'

class AbstractDictionary(object):
    """
    A base class for dictionary generation.

    Attributes:
        basis_func_number (int): The number of basis functions.
        B (numpy.ndarray): The B matrix.
    """

    def generate_B(self, inputs):
        """
        Generates the B matrix based on the inputs.

        Args:
            inputs (numpy.ndarray): The input array.

        Returns:
            numpy.ndarray: The generated B matrix.
        """
        target_dim = int(np.prod(inputs.shape[1:]))
        self.basis_func_number = self.n_psi_train + target_dim + 1
        # Form B matrix for flattened state coordinates.
        self.B = np.zeros((self.basis_func_number, target_dim))
        for i in range(0, target_dim):
            self.B[i + 1][i] = 1
        return self.B

class KoopmanNN(tf.keras.layers.Layer):
    def __init__(self, layer_sizes=[64, 64], n_psi_train=22, **kwargs):
        super(KoopmanNN, self).__init__(**kwargs)
        self.layer_sizes = layer_sizes
        self.n_psi_train = n_psi_train  # Using n_psi_train directly, consistent with DicNN
        self.input_layer = tf.keras.layers.Dense(layer_sizes[0], use_bias=False)
        self.hidden_layers = [tf.keras.layers.Dense(size, activation='tanh') for size in layer_sizes]
        self.output_layer = tf.keras.layers.Dense(n_psi_train)

    def call(self, inputs):
        raw_flat = tf.reshape(inputs, [tf.shape(inputs)[0], -1])
        x = self.input_layer(raw_flat)
        for layer in self.hidden_layers:
            x = layer(x)
        psi_x_train = self.output_layer(x)

        # Keep the constant term and raw coordinates as in the original ResKoopNet design.
        constant = tf.ones((tf.shape(inputs)[0], 1), dtype=inputs.dtype)
        outputs = tf.concat([constant, raw_flat, psi_x_train], axis=-1)
        return outputs

    def get_config(self):
        config = super(KoopmanNN, self).get_config()
        config.update({
            'layer_sizes': self.layer_sizes,
            'n_psi_train': self.n_psi_train
        })
        return config
    
    def generate_B(self, inputs):
        """
        Correctly generates the B matrix based on the inputs, using the proper attribute.
        """
        target_dim = int(np.prod(inputs.shape[1:]))
        # Use the flattened state dimension for non-vector inputs.
        self.basis_func_number = self.n_psi_train + target_dim + 1
        self.B = np.zeros((self.basis_func_number, target_dim))
        for i in range(0, target_dim):
            self.B[i + 1][i] = 1
        return self.B


class KoopmanCNN(tf.keras.layers.Layer):
    def __init__(
            self,
            input_shape,
            conv_filters=(32, 64, 128),
            kernel_size=3,
            pool_size=2,
            dense_units=(128,),
            n_psi_train=22,
            include_raw_state=True,
            **kwargs):
        super(KoopmanCNN, self).__init__(**kwargs)
        self.input_shape_spec = tuple(input_shape)
        self.conv_filters = tuple(conv_filters)
        self.kernel_size = kernel_size
        self.pool_size = pool_size
        self.dense_units = tuple(dense_units)
        self.n_psi_train = n_psi_train
        self.include_raw_state = include_raw_state

        # The last axis is assumed to be channel.
        self.spatial_ndim = len(self.input_shape_spec) - 1
        if self.spatial_ndim not in (1, 2, 3):
            raise ValueError(
                "KoopmanCNN expects input_shape to include a channel axis and be "
                "either (L, C), (H, W, C), or (D, H, W, C)."
            )

        if self.spatial_ndim == 1:
            conv_cls = Conv1D
            pool_cls = MaxPool1D
            gap_cls = GlobalAveragePooling1D
        elif self.spatial_ndim == 2:
            conv_cls = Conv2D
            pool_cls = MaxPool2D
            gap_cls = GlobalAveragePooling2D
        else:
            conv_cls = Conv3D
            pool_cls = MaxPool3D
            gap_cls = GlobalAveragePooling3D

        self.raw_flatten = Flatten(name='raw_flatten')
        self.conv_blocks = []
        for i, filters in enumerate(self.conv_filters):
            block = tf.keras.Sequential([
                conv_cls(
                    filters=filters,
                    kernel_size=self.kernel_size,
                    padding='same',
                    use_bias=False,
                    name=f'conv_{i}'
                ),
                BatchNormalization(name=f'bn_{i}'),
                Activation('relu', name=f'relu_{i}'),
                pool_cls(pool_size=self.pool_size, name=f'pool_{i}')
            ], name=f'conv_block_{i}')
            self.conv_blocks.append(block)

        self.global_pool = gap_cls(name='global_pool')
        self.hidden_layers = [
            Dense(size, activation='tanh', name=f'head_dense_{i}')
            for i, size in enumerate(self.dense_units)
        ]
        self.output_layer = Dense(n_psi_train, name='psi_train')

    def call(self, inputs):
        x = inputs
        for block in self.conv_blocks:
            x = block(x)
        x = self.global_pool(x)
        for layer in self.hidden_layers:
            x = layer(x)
        psi_x_train = self.output_layer(x)

        constant = tf.ones((tf.shape(inputs)[0], 1), dtype=inputs.dtype)
        if self.include_raw_state:
            raw_flat = self.raw_flatten(inputs)
            outputs = tf.concat([constant, raw_flat, psi_x_train], axis=-1)
        else:
            outputs = tf.concat([constant, psi_x_train], axis=-1)
        return outputs

    def get_config(self):
        config = super(KoopmanCNN, self).get_config()
        config.update({
            'input_shape': self.input_shape_spec,
            'conv_filters': self.conv_filters,
            'kernel_size': self.kernel_size,
            'pool_size': self.pool_size,
            'dense_units': self.dense_units,
            'n_psi_train': self.n_psi_train,
            'include_raw_state': self.include_raw_state
        })
        return config

    def generate_B(self, inputs):
        if not self.include_raw_state:
            raise ValueError(
                "compute_mode() requires explicit state coordinates. "
                "Set include_raw_state=True, or rewrite compute_mode to use a "
                "low-dimensional linear basis / decoder."
            )

        target_dim = int(np.prod(inputs.shape[1:]))
        self.basis_func_number = 1 + target_dim + self.n_psi_train
        self.B = np.zeros((self.basis_func_number, target_dim))
        self.B[1:1 + target_dim, :] = np.eye(target_dim)
        return self.B


class KoopmanAuxCNN(tf.keras.layers.Layer):
    def __init__(
            self,
            spec_input_shape,
            aux_dim,
            conv_filters=(32, 64, 128),
            kernel_size=3,
            pool_size=2,
            dense_units=(128,),
            n_psi_train=22,
            **kwargs):
        super(KoopmanAuxCNN, self).__init__(**kwargs)
        self.spec_input_shape = tuple(spec_input_shape)
        self.aux_dim = int(aux_dim)
        self.conv_filters = tuple(conv_filters)
        self.kernel_size = kernel_size
        self.pool_size = pool_size
        self.dense_units = tuple(dense_units)
        self.n_psi_train = n_psi_train

        self.spatial_ndim = len(self.spec_input_shape) - 1
        if self.spatial_ndim not in (1, 2, 3):
            raise ValueError(
                "KoopmanAuxCNN expects spec_input_shape to include a channel axis and be "
                "either (L, C), (H, W, C), or (D, H, W, C)."
            )

        if self.spatial_ndim == 1:
            conv_cls = Conv1D
            pool_cls = MaxPool1D
            gap_cls = GlobalAveragePooling1D
        elif self.spatial_ndim == 2:
            conv_cls = Conv2D
            pool_cls = MaxPool2D
            gap_cls = GlobalAveragePooling2D
        else:
            conv_cls = Conv3D
            pool_cls = MaxPool3D
            gap_cls = GlobalAveragePooling3D

        self.conv_blocks = []
        for i, filters in enumerate(self.conv_filters):
            block = tf.keras.Sequential([
                conv_cls(
                    filters=filters,
                    kernel_size=self.kernel_size,
                    padding='same',
                    use_bias=False,
                    name=f'conv_{i}'
                ),
                BatchNormalization(name=f'bn_{i}'),
                Activation('relu', name=f'relu_{i}'),
                pool_cls(pool_size=self.pool_size, name=f'pool_{i}')
            ], name=f'conv_block_{i}')
            self.conv_blocks.append(block)

        self.global_pool = gap_cls(name='global_pool')
        self.hidden_layers = [
            Dense(size, activation='tanh', name=f'head_dense_{i}')
            for i, size in enumerate(self.dense_units)
        ]
        self.output_layer = Dense(n_psi_train, name='psi_train')

    def _split_inputs(self, inputs):
        if isinstance(inputs, dict):
            if 'spec' not in inputs or 'aux' not in inputs:
                raise KeyError("KoopmanAuxCNN dict input must contain 'spec' and 'aux'.")
            spec = inputs['spec']
            aux = inputs['aux']
        elif isinstance(inputs, (list, tuple)) and len(inputs) == 2:
            spec, aux = inputs
        else:
            raise TypeError(
                "KoopmanAuxCNN expects inputs to be either a dict with keys "
                "'spec'/'aux' or a two-item tuple/list."
            )
        return spec, aux

    def call(self, inputs):
        spec, aux = self._split_inputs(inputs)

        x = spec
        for block in self.conv_blocks:
            x = block(x)
        x = self.global_pool(x)
        for layer in self.hidden_layers:
            x = layer(x)
        psi_x_train = self.output_layer(x)

        constant = tf.ones((tf.shape(spec)[0], 1), dtype=spec.dtype)
        aux = tf.cast(aux, spec.dtype)
        outputs = tf.concat([constant, aux, psi_x_train], axis=-1)
        return outputs

    def get_config(self):
        config = super(KoopmanAuxCNN, self).get_config()
        config.update({
            'spec_input_shape': self.spec_input_shape,
            'aux_dim': self.aux_dim,
            'conv_filters': self.conv_filters,
            'kernel_size': self.kernel_size,
            'pool_size': self.pool_size,
            'dense_units': self.dense_units,
            'n_psi_train': self.n_psi_train
        })
        return config

    def generate_B(self, inputs):
        _, aux = self._split_inputs(inputs)
        aux_shape = np.shape(aux)
        if len(aux_shape) < 2:
            raise ValueError("Auxiliary state must be shaped [batch, aux_dim].")

        target_dim = int(aux_shape[-1])
        if target_dim != self.aux_dim:
            raise ValueError(
                f"Auxiliary state dim ({target_dim}) does not match aux_dim ({self.aux_dim})."
            )

        self.basis_func_number = 1 + target_dim + self.n_psi_train
        self.B = np.zeros((self.basis_func_number, target_dim))
        self.B[1:1 + target_dim, :] = np.eye(target_dim)
        return self.B


class KoopmanSolver(object):


    '''
    Build the Koopman solver

    This part represents a Koopman solver that can be used to build and solve Koopman operator models.

    Attributes:
        dic (class): The dictionary class used for Koopman operator approximation.
        dic_func (function): The dictionary functions used for Koopman operator approximation.
        target_dim (int): The dimension of the variable of the equation.
        reg (float, optional): The regularization parameter when computing K. Defaults to 0.0.
        psi_x (None): Placeholder for the feature matrix of the input data.
        psi_y (None): Placeholder for the feature matrix of the output data.
    '''

    def __init__(
            self,
            dic,
            target_dim,
            reg=0.01,
            training_policy=None,
            analysis_dtype='float32',
            gram_dtype='float64',
            spectral_dtype='float64',
            residual_form='projected_kv'):
        """Initializer

        :param dic: dictionary
        :type dic: class
        :param target_dim: dimension of the variable of the equation
        :type target_dim: int
        :param reg: the regularization parameter when computing K, defaults to 0.0
        :type reg: float, optional
        """
        self.dic = dic  # dictionary class
        self.dic_func = dic.call  # dictionary functions
        if isinstance(target_dim, (tuple, list, tf.TensorShape, np.ndarray)):
            self.input_shape = tuple(int(dim) for dim in target_dim)
        else:
            self.input_shape = (int(target_dim),)
        self.target_dim = int(np.prod(self.input_shape))
        self.reg = reg
        self.residual_form = residual_form
        if self.residual_form not in ('projected_kv', 'projected_vlambda'):
            raise ValueError(
                "residual_form must be 'projected_kv' or 'projected_vlambda'."
            )
        self.training_policy = training_policy or (
            'mixed_float16' if tf.config.list_physical_devices('GPU') else 'float32'
        )
        self.analysis_dtype = tf.as_dtype(analysis_dtype)
        self.gram_dtype = tf.as_dtype(gram_dtype)
        self.spectral_dtype = np.dtype(spectral_dtype)
        self.training_input_dtype = tf.float32
        self.training_real_dtype = tf.float32
        self.training_complex_dtype = tf.complex64
        self.psi_x = None
        self.psi_y = None
        self.model = None
        self.layer_k = None
        self.analysis_dic = None
        self.output_dim = None
        self.eigenvectors_var = None
        self.eigenvalues_var = None
        self.eigenvectors_inv_B = None
        self.Psi_X = None
        self.Psi_Y = None
        self.outer_history = []
        self.eigvec_cond = None
        self.best_val_metric = np.nan
        self.best_outer_epoch = None
        self.best_outer_index = None
        self.best_train_metric = np.nan
        self.best_eigvec_cond = np.nan
        self.best_lr = np.nan
        self.best_reg = np.nan
        self.best_outer_summary = {}
        self.best_checkpoint_saved = False
        self.best_checkpoint_path = None
        self.final_checkpoint_path = None
        self.default_checkpoint_root = DEFAULT_CHECKPOINT_ROOT
        self.input_structure = None
        mixed_precision.set_global_policy(self.training_policy)

    def separate_data(self, data):
        data_x = data[0]
        data_y = data[1]
        return data_x, data_y

    def _is_nested_data(self, data):
        return isinstance(data, (dict, list, tuple))

    def _flatten_data(self, data):
        if self._is_nested_data(data):
            return tf.nest.flatten(data)
        return [data]

    def _data_length(self, data):
        first_leaf = self._flatten_data(data)[0]
        return int(np.shape(first_leaf)[0])

    def _slice_data(self, data, start=None, end=None, indices=None):
        def _slice_leaf(leaf):
            if indices is not None:
                return leaf[indices]
            return leaf[start:end]

        if self._is_nested_data(data):
            return tf.nest.map_structure(_slice_leaf, data)
        return _slice_leaf(data)

    def _convert_data_to_tensor(self, data, dtype):
        if self._is_nested_data(data):
            return tf.nest.map_structure(
                lambda leaf: tf.convert_to_tensor(leaf, dtype=dtype),
                data
            )
        return tf.convert_to_tensor(data, dtype=dtype)

    def _describe_data_shape(self, data):
        if self._is_nested_data(data):
            return tf.nest.map_structure(
                lambda leaf: tuple(int(dim) for dim in np.shape(leaf)),
                data
            )
        return tuple(int(dim) for dim in np.shape(data))

    def _infer_input_structure(self, data):
        if self._is_nested_data(data):
            return tf.nest.map_structure(
                lambda leaf: tuple(int(dim) for dim in np.shape(leaf)[1:]),
                data
            )
        return tuple(int(dim) for dim in np.shape(data)[1:])

    def _matches_unbatched_structure(self, data):
        if not self._is_nested_data(data) or self.input_structure is None:
            return False

        leaf_shapes = [
            tuple(int(dim) for dim in np.shape(leaf))
            for leaf in self._flatten_data(data)
        ]
        expected_shapes = list(tf.nest.flatten(self.input_structure))

        if len(leaf_shapes) != len(expected_shapes):
            return False
        return all(shape == expected for shape, expected in zip(leaf_shapes, expected_shapes))

    def _expand_batch_dim(self, data):
        if self._is_nested_data(data):
            return tf.nest.map_structure(
                lambda leaf: np.expand_dims(np.asarray(leaf), axis=0),
                data
            )
        return np.expand_dims(np.asarray(data), axis=0)

    def _batch_size_from_data(self, data):
        first_leaf = self._flatten_data(data)[0]
        return tf.shape(first_leaf)[0]

    def _build_input_placeholders(self, sample, prefix):
        if isinstance(sample, dict):
            return {
                key: self._build_input_placeholders(value, f"{prefix}_{key}")
                for key, value in sample.items()
            }
        if isinstance(sample, list):
            return [
                self._build_input_placeholders(value, f"{prefix}_{idx}")
                for idx, value in enumerate(sample)
            ]
        if isinstance(sample, tuple):
            return tuple(
                self._build_input_placeholders(value, f"{prefix}_{idx}")
                for idx, value in enumerate(sample)
            )

        input_shape = tuple(int(dim) for dim in np.shape(sample)[1:])
        return Input(input_shape, dtype=self.training_input_dtype, name=prefix)

    def _normalize_writable_path(self, path_value, default_root):
        """Normalize writable paths so relative paths land on the AutoDL data disk."""
        if path_value is None:
            return None

        expanded_path = os.path.expanduser(path_value)
        if os.path.isabs(expanded_path):
            normalized_path = os.path.normpath(expanded_path)
        else:
            normalized_path = os.path.normpath(os.path.join(default_root, expanded_path))

        read_only_root = os.path.normpath(READ_ONLY_AUTODL_ROOT)
        if normalized_path == read_only_root or normalized_path.startswith(read_only_root + os.sep):
            raise ValueError("Writable paths under /root/autodl-pub are not allowed.")

        return normalized_path

    def compute_final_info(self, reg_final=None):
        # Compute K
        reg_value = self.reg if reg_final is None else reg_final
        self.K = self.compute_K(self.dic_func,
                                self.data_x_train,
                                self.data_y_train,
                                reg=reg_value)
        self.eig_decomp(self.K)
        self.compute_mode()

    def eig_decomp(self, K):
        """ eigen-decomp of K """
        eig_start = time.time()
        print(f"eig_decomp: start, K shape={np.shape(K)}, spectral_dtype={self.spectral_dtype}")
        tensor_to_numpy_start = time.time()
        K_np = np.array(K.numpy(), dtype=self.spectral_dtype, order='F', copy=True)
        print(
            "eig_decomp: tensor->numpy done in "
            f"{time.time() - tensor_to_numpy_start:.2f}s, "
            f"norm={np.linalg.norm(K_np):.3e}, "
            f"max_abs={np.max(np.abs(K_np)):.3e}, "
            f"has_nan={np.isnan(K_np).any()}, has_inf={np.isinf(K_np).any()}"
        )
        eigsolve_start = time.time()
        print("eig_decomp: starting scipy.linalg.eig")
        self.eigenvalues, self.eigenvectors = scipy.linalg.eig(
            K_np,
            left=False,
            right=True,
            overwrite_a=True,
            check_finite=False
        )
        eigsolve_elapsed = time.time() - eigsolve_start
        print(f"eig_decomp: scipy.linalg.eig done in {eigsolve_elapsed:.2f}s")
        idx = self.eigenvalues.real.argsort()[::-1]
        self.eigenvalues = self.eigenvalues[idx]
        self.eigenvectors = self.eigenvectors[:, idx]
        self.eigvec_cond = float(np.linalg.cond(self.eigenvectors))
        print(f"eig_decomp: eigenvector condition number={self.eigvec_cond:.3e}")
        print(f"eig_decomp: done in {time.time() - eig_start:.2f}s")

    def _ensure_batched_input(self, data_x):
        if self._is_nested_data(data_x):
            if self._matches_unbatched_structure(data_x):
                data_x = self._expand_batch_dim(data_x)
            return data_x
        data_x = np.asarray(data_x)
        if tuple(data_x.shape) == self.input_shape:
            data_x = np.expand_dims(data_x, axis=0)
        return data_x

    def eigenfunctions(self, data_x):
        """ estimated eigenfunctions """
        data_x = self._ensure_batched_input(data_x)
        analysis_dic = self._refresh_analysis_dictionary(self._slice_data(data_x, 0, 1))
        psi_x = analysis_dic(self._convert_data_to_tensor(data_x, self.analysis_dtype)).numpy()
        psi_x = psi_x.astype(self.spectral_dtype, copy=False)
        val = np.matmul(psi_x, self.eigenvectors)
        return val

    def compute_mode(self):
        self.basis_func_number = self.K.shape[0]

        # Form B matrix
        self.B = self.dic.generate_B(self.data_x_train)

        # Compute modes
        solve_start = time.time()
        print("compute_mode: starting solve(V, B)")
        self.eigenvectors_inv_B = np.linalg.solve(self.eigenvectors, self.B)
        print(f"compute_mode: solve(V, B) done in {time.time() - solve_start:.2f}s")
        self.modes = self.eigenvectors_inv_B.T
        return self.modes

    def calc_psi_next(self, data_x, K):
        psi_x = self.dic_func(data_x)
        psi_next = tf.matmul(psi_x, K)
        return psi_next

    def predict(self, x0, traj_len):
        """ predict the trajectory """
        if self._is_nested_data(x0):
            raise NotImplementedError(
                "predict() is not supported for multimodal inputs because the "
                "current solver only reconstructs the explicit low-dimensional state."
            )
        x_curr = self._ensure_batched_input(x0)
        traj = [x_curr]
        for _ in range(traj_len - 1):
            efunc = self.eigenfunctions(x_curr)
            x_next = np.matmul(self.modes, (self.eigenvalues * efunc).T)
            x_curr = (x_next.real).T.reshape((x_next.shape[1],) + self.input_shape)
            traj.append(x_curr)
        traj = np.stack(traj, axis=1)
        return traj.squeeze()

    def _refresh_analysis_dictionary(self, sample):
        training_sample = self._convert_data_to_tensor(sample, self.training_input_dtype)
        _ = self.dic(training_sample)

        dic_cls = self.dic.__class__
        if not hasattr(self.dic, 'get_config'):
            raise TypeError("Dictionary must implement get_config() to build the analysis copy.")

        current_policy = mixed_precision.global_policy()
        mixed_precision.set_global_policy(self.analysis_dtype.name)
        try:
            if hasattr(dic_cls, 'from_config'):
                analysis_dic = dic_cls.from_config(self.dic.get_config())
            else:
                analysis_dic = dic_cls(**self.dic.get_config())

            analysis_sample = self._convert_data_to_tensor(sample, self.analysis_dtype)
            _ = analysis_dic(analysis_sample)
        finally:
            mixed_precision.set_global_policy(current_policy.name)

        cast_weights = [
            weight.astype(self.analysis_dtype.as_numpy_dtype, copy=False)
            for weight in self.dic.get_weights()
        ]
        analysis_dic.set_weights(cast_weights)
        self.analysis_dic = analysis_dic
        return self.analysis_dic

    # def compute_K(self, dic, data_x, data_y, reg):
    #     psi_x = dic(data_x)
    #     psi_y = dic(data_y)
    #     # Compute Psi_X and Psi_Y
    #     self.Psi_X = dic(data_x)
    #     self.Psi_Y = dic(data_y)
    #     psi_xt = tf.transpose(psi_x)
    #     idmat = tf.eye(psi_x.shape[-1], dtype='float64')
    #     xtx_inv = tf.linalg.pinv(reg * idmat + tf.matmul(psi_xt, psi_x))
    #     xty = tf.matmul(psi_xt, psi_y)
    #     self.K_reg = tf.matmul(xtx_inv, xty)
    #     return self.K_reg
    # def compute_K(self, dic, data_x, data_y, reg):
    #     print("compute_K: start")
    #     print(f"compute_K: data_x shape={data_x.shape}, data_y shape={data_y.shape}")
    #     # build dictionary matrices
    #     psi_x = self.batch_apply_dic(data_x, self.batch_size)
    #     print(f"compute_K: psi_x computed, shape={psi_x.shape}")
    #     psi_y = self.batch_apply_dic(data_y, self.batch_size)
    #     print(f"compute_K: psi_y computed, shape={psi_y.shape}")

    #     # store for later inspection
    #     print("compute_K: computing and storing Psi_X, Psi_Y")
    #     self.Psi_X = psi_x 
    #     self.Psi_Y = psi_y
    #     print(f"compute_K: Psi_X shape={self.Psi_X.shape}, Psi_Y shape={self.Psi_Y.shape}")

    #     # form normal equations
    #     print("compute_K: forming normal equation components")
    #     psi_xt = tf.transpose(psi_x)
    #     print(psi_xt.shape, psi_x.shape)
    #     idmat = tf.eye(psi_x.shape[-1], dtype='float64')

    #     print("compute_K: computing (X^T X + reg*I)^{-1}")
    #     xtx_inv = tf.linalg.pinv(reg * idmat + tf.matmul(psi_xt, psi_x))
    #     print("compute_K: xtx_inv computed")

    #     print("compute_K: computing X^T Y")
    #     xty = tf.matmul(psi_xt, psi_y)
    #     print(f"compute_K: xty computed, shape={xty.shape}")

    #     # compute and return Koopman matrix
    #     print("compute_K: computing K_reg")
    #     self.K_reg = tf.matmul(xtx_inv, xty)
    #     print(f"compute_K: done, K_reg shape={self.K_reg.shape}")
    #     return self.K_reg

    # def batch_apply_dic(self, data, batch_size=128):
    #     parts = []
    #     for i in range(0, data.shape[0], batch_size):
    #         print(i)
    #         parts.append(self.dic_func(data[i:i+batch_size]))
    #     return tf.concat(parts, axis=0)

    # def compute_K(self, dic, data_x, data_y, reg):
    #     print("compute_K: start (streaming Gram)")
    #     d = dic(data_x[:1]).shape[-1]   # dictionary dim
    #     XtX = tf.zeros((d, d), dtype='float64')
    #     XtY = tf.zeros((d, d), dtype='float64')

    #     for i in range(0, data_x.shape[0], self.batch_size):
    #         Xb = data_x[i:i+self.batch_size]
    #         Yb = data_y[i:i+self.batch_size]

    #         psi_x_b = dic(Xb)   # [b, d]
    #         psi_y_b = dic(Yb)   # [b, d]

    #         # accumulate
    #         XtX += tf.matmul(psi_x_b, psi_x_b, transpose_a=True)
    #         XtY += tf.matmul(psi_x_b, psi_y_b, transpose_a=True)

    #         print(f"  batch {i} → XtX trace={tf.linalg.trace(XtX):.2e}")

    #     # now form and solve normal equations
    #     idmat = tf.eye(d, dtype='float64')
    #     inv_term = tf.linalg.pinv(XtX + reg * idmat)
    #     K_reg = tf.matmul(inv_term, XtY)

    #     print(f"compute_K: done, K_reg shape={K_reg.shape}")
    #     return K_reg


    def compute_K(self, dic, data_x, data_y, reg):
        """
        Streaming computation of the Koopman matrix with a progress bar.
        """
        compute_k_start = time.time()
        sample = self._slice_data(data_x, 0, 1)
        analysis_dic = self._refresh_analysis_dictionary(sample)
        d = analysis_dic(self._convert_data_to_tensor(sample, self.analysis_dtype)).shape[-1]

        # initialize accumulators
        XtX_sum = tf.zeros((d, d), dtype=self.gram_dtype)
        XtY_sum = tf.zeros((d, d), dtype=self.gram_dtype)
        n_total = 0

        total = self._data_length(data_x)

        print(self.batch_size, total)
        # iterate in batches with tqdm progress bar
        for start in tqdm(range(0, total, self.batch_size), desc="compute_K batches", unit="batch"):
            end = min(start + self.batch_size, total)
            Xb = self._slice_data(data_x, start, end)
            Yb = self._slice_data(data_y, start, end)

            psi_x_b = analysis_dic(self._convert_data_to_tensor(Xb, self.analysis_dtype))
            psi_y_b = analysis_dic(self._convert_data_to_tensor(Yb, self.analysis_dtype))
            psi_x_b = tf.cast(psi_x_b, self.gram_dtype)
            psi_y_b = tf.cast(psi_y_b, self.gram_dtype)
            batch_n = int(end - start)

            # accumulate Gram and cross-covariance
            XtX_sum += tf.matmul(psi_x_b, psi_x_b, transpose_a=True)
            XtY_sum += tf.matmul(psi_x_b, psi_y_b, transpose_a=True)
            n_total += batch_n

        # form and solve normal equations
        if n_total == 0:
            raise ValueError("compute_K received an empty dataset.")

        n_total_tensor = tf.cast(n_total, self.gram_dtype)
        G = XtX_sum / n_total_tensor
        A = XtY_sum / n_total_tensor
        xtx_diag = tf.linalg.diag_part(G)
        xtx_diag_min = float(tf.reduce_min(xtx_diag).numpy())
        xtx_diag_max = float(tf.reduce_max(xtx_diag).numpy())
        xtx_has_nan = bool(tf.reduce_any(tf.math.is_nan(G)).numpy())
        xtx_has_inf = bool(tf.reduce_any(tf.math.is_inf(G)).numpy())
        print(
            "compute_K: batches done, "
            f"d={d}, n_total={n_total}, G diag min={xtx_diag_min:.3e}, G diag max={xtx_diag_max:.3e}, "
            f"has_nan={xtx_has_nan}, has_inf={xtx_has_inf}, elapsed={time.time() - compute_k_start:.2f}s"
        )
        idmat = tf.eye(d, dtype=self.gram_dtype)
        system = G + tf.cast(reg, self.gram_dtype) * idmat
        solve_start = time.time()
        print("compute_K: starting linear solve")
        use_pinv = False
        try:
            K_reg = tf.linalg.solve(system, A)
            if not bool(tf.reduce_all(tf.math.is_finite(K_reg)).numpy()):
                use_pinv = True
                print("compute_K: linear solve produced non-finite values, falling back to pinv")
        except (tf.errors.InvalidArgumentError, tf.errors.InternalError) as exc:
            use_pinv = True
            print(f"compute_K: linear solve failed with {type(exc).__name__}, falling back to pinv")

        if use_pinv:
            pinv_start = time.time()
            print("compute_K: starting pinv fallback")
            inv_term = tf.linalg.pinv(system)
            print(f"compute_K: pinv fallback done in {time.time() - pinv_start:.2f}s")
            matmul_start = time.time()
            print("compute_K: starting final matmul")
            K_reg = tf.matmul(inv_term, A)
            print(
                f"compute_K: final matmul done in {time.time() - matmul_start:.2f}s, "
                f"total compute_K time={time.time() - compute_k_start:.2f}s"
            )
        else:
            print(
                f"compute_K: linear solve done in {time.time() - solve_start:.2f}s, "
                f"total compute_K time={time.time() - compute_k_start:.2f}s"
            )

        return K_reg




    def get_Psi_X(self):
        if self.Psi_X is None:
            if not hasattr(self, 'data_x_train'):
                raise ValueError("Training data is not available for Psi_X.")
            self.Psi_X = self._compute_dictionary_matrix(self.data_x_train)
        return self.Psi_X

    def get_Psi_Y(self):
        if self.Psi_Y is None:
            if not hasattr(self, 'data_y_train'):
                raise ValueError("Training data is not available for Psi_Y.")
            self.Psi_Y = self._compute_dictionary_matrix(self.data_y_train)
        return self.Psi_Y

    def _compute_dictionary_matrix(self, data):
        if self._data_length(data) == 0:
            return np.empty((0, 0), dtype=self.spectral_dtype)

        analysis_dic = self._refresh_analysis_dictionary(self._slice_data(data, 0, 1))
        psi = analysis_dic(self._convert_data_to_tensor(data, self.analysis_dtype)).numpy()
        return psi.astype(self.spectral_dtype, copy=False)

    def _pack_complex_residual(self, residual):
        return tf.concat(
            [tf.math.real(residual), tf.math.imag(residual)],
            axis=-1
        )

    def _packed_residual_loss(self, y_true, y_pred):
        diff = tf.cast(y_pred, self.training_real_dtype) - tf.cast(y_true, self.training_real_dtype)
        return tf.reduce_mean(tf.reduce_sum(tf.square(diff), axis=-1))

    def _build_training_residual(self, psi_x, psi_y):
        if self.eigenvectors_var is None or self.eigenvalues_var is None:
            raise ValueError("Spectral variables must be initialized before build_model().")

        if self.residual_form == 'projected_kv':
            if self.layer_k is None:
                raise ValueError("Layer_K must be initialized for residual_form='projected_kv'.")
            psi_x_real = tf.cast(psi_x, self.training_real_dtype)
            psi_y_real = tf.cast(psi_y, self.training_real_dtype)
            psi_next = self.layer_k(psi_x_real)
            psi_diff = tf.cast(psi_y_real - psi_next, self.training_complex_dtype)
            return tf.matmul(psi_diff, self.eigenvectors_var)

        psi_x_complex = tf.cast(psi_x, self.training_complex_dtype)
        psi_y_complex = tf.cast(psi_y, self.training_complex_dtype)
        psi_y_v = tf.matmul(psi_y_complex, self.eigenvectors_var)
        psi_x_v = tf.matmul(psi_x_complex, self.eigenvectors_var)
        return psi_y_v - psi_x_v * self.eigenvalues_var

    def _sync_training_spectral_state(self):
        if self.eigenvectors_var is not None:
            self.eigenvectors_var.assign(tf.cast(self.eigenvectors, self.training_complex_dtype))
        if self.eigenvalues_var is not None:
            self.eigenvalues_var.assign(tf.cast(self.eigenvalues, self.training_complex_dtype))
        if self.residual_form == 'projected_kv' and self.model is not None:
            layer_k_weights = self.model.get_layer('Layer_K').weights[0]
            layer_k_weights.assign(tf.cast(self.K, layer_k_weights.dtype))

    def evaluate_spectral_residual(self, data_x, data_y, residual_form=None):
        residual_form = residual_form or self.residual_form
        total_sq = 0.0
        total_n = 0
        batch_size = getattr(self, 'batch_size', 1024)
        eigvecs = tf.cast(self.eigenvectors, self.training_complex_dtype)
        eigvals = tf.cast(self.eigenvalues, self.training_complex_dtype)
        k_matrix = tf.cast(self.K, self.training_real_dtype)

        total = self._data_length(data_x)
        for start in range(0, total, batch_size):
            end = min(start + batch_size, total)
            x_batch = self._convert_data_to_tensor(
                self._slice_data(data_x, start, end),
                self.training_input_dtype
            )
            y_batch = self._convert_data_to_tensor(
                self._slice_data(data_y, start, end),
                self.training_input_dtype
            )
            psi_x = self.dic_func(x_batch)
            psi_y = self.dic_func(y_batch)

            if residual_form == 'projected_kv':
                psi_next = tf.matmul(tf.cast(psi_x, self.training_real_dtype), k_matrix)
                residual = tf.matmul(
                    tf.cast(tf.cast(psi_y, self.training_real_dtype) - psi_next, self.training_complex_dtype),
                    eigvecs
                )
            else:
                psi_x_v = tf.matmul(tf.cast(psi_x, self.training_complex_dtype), eigvecs)
                psi_y_v = tf.matmul(tf.cast(psi_y, self.training_complex_dtype), eigvecs)
                residual = psi_y_v - psi_x_v * eigvals

            packed = self._pack_complex_residual(residual)
            batch_sq = tf.reduce_sum(tf.square(packed), axis=-1)
            total_sq += float(tf.reduce_sum(batch_sq).numpy())
            total_n += int(end - start)

        if total_n == 0:
            return 0.0
        return total_sq / total_n


    '''
    Build the Koopman model with dictionary learning
    '''
    def build_model(self):
        """Build model with trainable dictionary

        The loss function is ||Psi(y) - K Psi(x)||^2 .

        """
        if self._is_nested_data(self.data_x_train):
            inputs_x = self._build_input_placeholders(self.data_x_train, 'input_x')
            inputs_y = self._build_input_placeholders(self.data_y_train, 'input_y')
        else:
            inputs_x = Input(self.input_shape, dtype=self.training_input_dtype)
            inputs_y = Input(self.input_shape, dtype=self.training_input_dtype)

        self.psi_x = self.dic_func(inputs_x)
        self.psi_y = self.dic_func(inputs_y)

        self.layer_k = None
        if self.residual_form == 'projected_kv':
            self.layer_k = Dense(units=self.psi_y.shape[-1],
                                 use_bias=False,
                                 name='Layer_K',
                                 trainable=False,
                                 dtype='float32')

        residual = self._build_training_residual(self.psi_x, self.psi_y)
        outputs = self._pack_complex_residual(residual)
        model = Model(inputs=[inputs_x, inputs_y], outputs=outputs)
        return model

    def train_psi(self, model, epochs, callbacks):
        """Train the trainable part of the dictionary

        :param model: koopman model
        :type model: model
        :param epochs: the number of training epochs before computing K for each inner training epoch
        :type epochs: int
        :return: history
        :rtype: history callback object
        """
        # history = model.fit(
        #     x=self.data_train,
        #     y=self.zeros_data_y_train,
        #     epochs=epochs,
        #     callbacks = callbacks,
        #     validation_data=(
        #         self.data_valid,
        #         self.zeros_data_y_valid),
        #     batch_size=self.batch_size,
        #     verbose=1)
        history = model.fit(
        x=self.ds_train,
        validation_data=self.ds_valid,
        epochs=epochs,
        callbacks=callbacks,
        verbose=1
        )

        return history

    def get_basis(self, x, y):
        """Returns the dictionary(matrix) consisting of basis.

        :param x: array of snapshots
        :type x: numpy array
        :param y:array of snapshots
        :type y: numpy array
        """
        psi_x = self.dic_func(x)
        # Calculate column norms
        psi_x_column_norms = np.linalg.norm(psi_x, axis=0)
        # Handle the case where norm is zero
        psi_x_column_norms[psi_x_column_norms == 0] = 1
        psi_x_normalized = psi_x / psi_x_column_norms

        # Repeat the steps for psi_y
        psi_y = self.dic_func(y)
        # Calculate column norms
        psi_y_column_norms = np.linalg.norm(psi_y, axis=0)
        # Handle the case where norm is zero
        psi_y_column_norms[psi_y_column_norms == 0] = 1
        psi_y_normalized = psi_y / psi_y_column_norms

        return psi_x_normalized, psi_y_normalized

    # def get_derivatives(self, inputs):
    #     # Compute the first and second derivatives
    #     inputs = tf.convert_to_tensor(inputs)
    #     with tf.GradientTape() as tape2:
    #         with tf.GradientTape() as tape1:
    #             tape1.watch(inputs)
    #             tape2.watch(inputs)
    #             outputs = self.dic(inputs)
    #         first_derivatives = tape1.batch_jacobian(outputs, inputs)
    #     second_derivatives = tape2.batch_jacobian(first_derivatives, inputs)
    #     # print(outputs.shape, first_derivatives.shape, second_derivatives.shape)
    #     # return (outputs, first_derivatives, second_derivatives)
    #     return (first_derivatives, second_derivatives)
    
# inside build(), after you have self.data_x_train & self.data_y_train
    def _make_ds(self, x_array, y_array):
        x_tensor = self._convert_data_to_tensor(x_array, self.training_input_dtype)
        y_tensor = self._convert_data_to_tensor(y_array, self.training_input_dtype)
        ds = tf.data.Dataset.from_tensor_slices((x_tensor, y_tensor))
        ds = ds.batch(self.batch_size)

        def map_to_zero(x_batch, y_batch):
            batch_size = self._batch_size_from_data(x_batch)
            zeros = tf.zeros((batch_size, self.output_dim), dtype=tf.float32)
            return (x_batch, y_batch), zeros

        ds = ds.map(map_to_zero, num_parallel_calls=tf.data.AUTOTUNE)
        return ds.prefetch(tf.data.AUTOTUNE)


    def build(
            self,
            data_train,
            data_valid,
            epochs,
            batch_size,
            lr,
            log_interval,
            lr_decay_factor,
            Nepoch,
            end_condition=1e-5,
            checkpoint_path=None,
            best_checkpoint_path=None,
            save_best_only=True,
            resume=False
        ):
        print("build: start")
        self.data_train = data_train
        self.batch_size = batch_size
        self.outer_history = []
        self.best_val_metric = np.nan
        self.best_outer_epoch = None
        self.best_outer_index = None
        self.best_train_metric = np.nan
        self.best_eigvec_cond = np.nan
        self.best_lr = np.nan
        self.best_reg = np.nan
        self.best_outer_summary = {}
        self.best_checkpoint_saved = False
        checkpoint_path = self._normalize_writable_path(
            checkpoint_path,
            self.default_checkpoint_root
        )
        best_checkpoint_path = self._normalize_writable_path(
            best_checkpoint_path,
            self.default_checkpoint_root
        )
        self.best_checkpoint_path = best_checkpoint_path
        self.final_checkpoint_path = checkpoint_path
        
        self.data_x_train, self.data_y_train = self.separate_data(self.data_train)
        print(
            "build: separate_data done, "
            f"x_train shape={self._describe_data_shape(self.data_x_train)}, "
            f"y_train shape={self._describe_data_shape(self.data_y_train)}"
        )
        self.input_structure = self._infer_input_structure(self.data_x_train)
        if not self._is_nested_data(self.data_x_train):
            observed_input_shape = tuple(int(dim) for dim in self.data_x_train.shape[1:])
            if observed_input_shape != self.input_shape:
                print(f"build: overriding input_shape from {self.input_shape} to {observed_input_shape}")
                self.input_shape = observed_input_shape
                self.target_dim = int(np.prod(self.input_shape))
        self.data_valid = data_valid

        print("build: compute_final_info start")
        self.compute_final_info()
        print(f"build: compute_final_info done, K shape={self.K.shape}")
        self.output_dim = int(self.K.shape[0]) * 2
        self.eigenvectors_var = tf.Variable(
            tf.cast(self.eigenvectors, self.training_complex_dtype),
            trainable=False,
            name='eigenvectors'
        )
        self.eigenvalues_var = tf.Variable(
            tf.cast(self.eigenvalues, self.training_complex_dtype),
            trainable=False,
            name='eigenvalues'
        )

        self.ds_train = self._make_ds(self.data_x_train, self.data_y_train)
        self.ds_valid = self._make_ds(self.data_valid[0], self.data_valid[1])

        # self.zeros_data_y_train = tf.zeros_like(self.dic_func(self.data_y_train))
        # self.zeros_data_y_valid = tf.zeros_like(self.dic_func(self.data_valid[1]))
        print(f"creating datasets")

        self.model = self.build_model()
        print("build: build_model done")

        opt = Adam(lr)
        self.model.compile(optimizer=opt, loss=self._packed_residual_loss)
        self._sync_training_spectral_state()

        K_var = tf.Variable(self.K, trainable=False, name='K')
        reg_var = tf.Variable(self.reg, trainable=False, name='reg')

        checkpoint = None
        latest_checkpoint = None
        checkpoint_dir = None
        best_checkpoint_dir = None
        if checkpoint_path or best_checkpoint_path:
            if checkpoint_path:
                checkpoint_dir = os.path.dirname(checkpoint_path) or '.'
            if best_checkpoint_path:
                best_checkpoint_dir = os.path.dirname(best_checkpoint_path) or '.'
            checkpoint = tf.train.Checkpoint(
                model=self.model,
                K=K_var,
                eigenvectors=self.eigenvectors_var,
                eigenvalues=self.eigenvalues_var,
                reg=reg_var
            )
        if checkpoint_path:
            if resume:
                latest_checkpoint = tf.train.latest_checkpoint(checkpoint_dir)
                if latest_checkpoint:
                    checkpoint.restore(latest_checkpoint)
                    self.K = K_var.numpy()
                    self.eigenvectors = self.eigenvectors_var.numpy().astype(np.complex128, copy=False)
                    self.eigenvalues = self.eigenvalues_var.numpy().astype(np.complex128, copy=False)
                    self.reg = float(reg_var.numpy())
                    self._sync_training_spectral_state()
                    print(f"build: restored checkpoint from {latest_checkpoint}")
                else:
                    print("build: no checkpoint found, starting fresh")
        if checkpoint_path and not resume:
            print("build: resume disabled, skipping checkpoint restore")
        print(f"build: optimizer lr={opt.lr.numpy()}")

        losses = []
        val_losses = []
        learning_rate_changes = []
        stop_flag = 0
        best_val_metric = np.inf
        best_outer_epoch = None
        best_checkpoint_saved = False

        def _sync_checkpoint_state():
            K_var.assign(tf.cast(self.K, K_var.dtype))
            self.eigenvectors_var.assign(tf.cast(self.eigenvectors, self.training_complex_dtype))
            self.eigenvalues_var.assign(tf.cast(self.eigenvalues, self.training_complex_dtype))
            reg_var.assign(tf.cast(self.reg, reg_var.dtype))

        for i in range(epochs):
            outer_epoch_start = time.time()
            print(f"build: Outer Epoch {i+1}/{epochs}")
            print("build: compute_K start")
            self.K = self.compute_K(self.dic_func, self.data_x_train, self.data_y_train, self.reg)
            self.eig_decomp(self.K)
            self.compute_mode()
            self._sync_training_spectral_state()
            print("build: compute_K done")

            if self.residual_form == 'projected_kv':
                print("build: assigned K to Layer_K")
            else:
                print("build: synced V/Lambda into training graph")

            print("build: train_psi start")
            self.history = self.train_psi(
                self.model,
                epochs=Nepoch,
                callbacks=[EarlyStopping(monitor='val_loss', min_delta=1e-6, patience=2, verbose=1, mode='min')]
            )
            print(f"build: train_psi done, history length={len(self.history.history['loss'])}")

            gc.collect()
            mem = memory_usage()[0]
            print(f"build: post-gc memory usage={mem} MiB")

            print("build: post-train compute_K start")
            self.K = self.compute_K(self.dic_func, self.data_x_train, self.data_y_train, self.reg)
            self.eig_decomp(self.K)
            self.compute_mode()
            self._sync_training_spectral_state()
            print("build: post-train compute_K done")

            inner_train_last = float(self.history.history['loss'][-1]) if self.history.history.get('loss') else np.nan
            inner_val_last = float(self.history.history['val_loss'][-1]) if self.history.history.get('val_loss') else np.nan
            train_metric = self.evaluate_spectral_residual(self.data_x_train, self.data_y_train)
            val_metric = self.evaluate_spectral_residual(self.data_valid[0], self.data_valid[1])
            self.outer_history.append({
                'outer_epoch': i + 1,
                'inner_train_last': inner_train_last,
                'inner_val_last': inner_val_last,
                'train_metric': train_metric,
                'val_metric': val_metric,
                'eigvec_cond': self.eigvec_cond,
                'lr': float(opt.lr.numpy()),
                'reg': float(self.reg)
            })
            print(
                f"build: outer metrics train={train_metric:.6e}, "
                f"val={val_metric:.6e}, eigvec_cond={self.eigvec_cond:.3e}"
            )

            current_outer_index = len(self.outer_history) - 1
            current_is_better = False
            if best_outer_epoch is None:
                current_is_better = True
            elif np.isfinite(val_metric) and (not np.isfinite(best_val_metric) or val_metric < best_val_metric):
                current_is_better = True

            if current_is_better:
                best_val_metric = float(val_metric)
                best_outer_epoch = i + 1
                self.best_val_metric = float(val_metric)
                self.best_outer_epoch = i + 1
                self.best_outer_index = current_outer_index
                self.best_train_metric = float(train_metric)
                self.best_eigvec_cond = float(self.eigvec_cond)
                self.best_lr = float(opt.lr.numpy())
                self.best_reg = float(self.reg)
                self.best_outer_summary = {
                    'best_val_metric': self.best_val_metric,
                    'best_outer_epoch': self.best_outer_epoch,
                    'best_outer_index': self.best_outer_index,
                    'best_train_metric': self.best_train_metric,
                    'best_eigvec_cond': self.best_eigvec_cond,
                    'best_lr': self.best_lr,
                    'best_reg': self.best_reg
                }
                print(
                    f"build: new best outer val_metric={self.best_val_metric:.6e} "
                    f"at epoch {self.best_outer_epoch}"
                )

            should_save_best_checkpoint = (
                checkpoint is not None
                and best_checkpoint_path is not None
                and (current_is_better or not save_best_only)
            )
            if should_save_best_checkpoint:
                if best_checkpoint_dir and not os.path.exists(best_checkpoint_dir):
                    os.makedirs(best_checkpoint_dir, exist_ok=True)
                _sync_checkpoint_state()
                checkpoint.save(file_prefix=best_checkpoint_path)
                best_checkpoint_saved = True
                self.best_checkpoint_saved = True
                if current_is_better:
                    print("build: best checkpoint saved")

            if i % log_interval == 0:
                losses.extend(self.history.history['loss'])
                val_losses.extend(self.history.history['val_loss'])

            if len(self.outer_history) > 1:
                prev_val_metric = self.outer_history[-2]['val_metric']
                if val_metric > prev_val_metric:
                    print("build: outer validation metric increased, decaying lr")
                    opt.lr.assign(opt.lr * lr_decay_factor)

            learning_rate_changes.append(opt.lr.numpy())
            print(f"build: Outer Epoch {i+1} finished in {time.time() - outer_epoch_start:.2f}s")

            if len(self.outer_history) >= 3:
                window = min(5, len(self.outer_history))
                recent_vals = [entry['val_metric'] for entry in self.outer_history[-window:]]
                avg_metric_change = np.mean(np.abs(np.diff(recent_vals)))
                if avg_metric_change < end_condition:
                    print("build: outer validation metric stabilized, stopping training")
                    stop_flag = 1
                    break

        print("build: final compute_final_info start")
        self.compute_final_info()
        print("build: final compute_final_info done")
        self._sync_training_spectral_state()
        _sync_checkpoint_state()
        self.best_checkpoint_saved = best_checkpoint_saved

        if checkpoint_path:
            if checkpoint_dir and not os.path.exists(checkpoint_dir):
                os.makedirs(checkpoint_dir, exist_ok=True)
            checkpoint.save(file_prefix=checkpoint_path)
            print("build: checkpoint saved")

        return losses, val_losses, stop_flag, learning_rate_changes
