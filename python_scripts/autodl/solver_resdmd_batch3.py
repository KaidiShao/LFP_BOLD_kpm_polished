import json
import os
# from autograd import jacobian, hessian
import time
import tensorflow as tf
from tensorflow.keras.layers import Dense, Layer, Concatenate, Input
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


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.generic,)):
        return value.item()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _save_json_file(path_value, payload):
    if not path_value:
        return

    os.makedirs(os.path.dirname(path_value), exist_ok=True)
    temp_path = f"{path_value}.tmp"
    with open(temp_path, "w", encoding="utf-8") as handle:
        json.dump(_json_ready(payload), handle, indent=2, sort_keys=True)
    os.replace(temp_path, path_value)


def _load_json_file(path_value):
    if not path_value or (not os.path.exists(path_value)):
        return None
    with open(path_value, "r", encoding="utf-8") as handle:
        return json.load(handle)

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
        target_dim = inputs.shape[-1]
        self.basis_func_number = self.n_psi_train + target_dim + 1
        # Form B matrix
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
        x = self.input_layer(inputs)
        for layer in self.hidden_layers:
            x = layer(x)
        psi_x_train = self.output_layer(x)
        
        # Directly integrating the generation of the constant and concatenation as done in PsiNN
        constant = tf.ones_like(inputs[:, :1])
        outputs = tf.keras.layers.Concatenate()([constant, inputs, psi_x_train])
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
        target_dim = inputs.shape[-1]
        # Use n_psi_train instead of n_dic_customized
        self.basis_func_number = self.n_psi_train + target_dim + 1
        self.B = np.zeros((self.basis_func_number, target_dim))
        for i in range(0, target_dim):
            self.B[i + 1][i] = 1
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
        self.target_dim = target_dim
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
        self.best_state_path = None
        self.final_state_path = None
        self.default_checkpoint_root = DEFAULT_CHECKPOINT_ROOT
        mixed_precision.set_global_policy(self.training_policy)

    def separate_data(self, data):
        data_x = data[0]
        data_y = data[1]
        return data_x, data_y

    def _normalize_writable_path(self, path_value, default_root):
        """Normalize writable paths so relative paths land on the local results disk."""
        if path_value is None:
            return None

        expanded_path = os.path.expanduser(path_value)
        if os.path.isabs(expanded_path):
            normalized_path = os.path.normpath(expanded_path)
        else:
            normalized_path = os.path.normpath(os.path.join(default_root, expanded_path))

        read_only_root = os.path.normpath(READ_ONLY_AUTODL_ROOT)
        if normalized_path == read_only_root or normalized_path.startswith(read_only_root + os.sep):
            raise ValueError(f"Writable paths under {READ_ONLY_AUTODL_ROOT} are not allowed.")

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

    def eigenfunctions(self, data_x):
        """ estimated eigenfunctions """
        analysis_dic = self._refresh_analysis_dictionary(data_x[:1])
        psi_x = analysis_dic(tf.convert_to_tensor(data_x, dtype=self.analysis_dtype)).numpy()
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
        traj = [x0]
        for _ in range(traj_len - 1):
            x_curr = traj[-1]
            efunc = self.eigenfunctions(x_curr)
            x_next = np.matmul(self.modes, (self.eigenvalues * efunc).T)
            traj.append((x_next.real).T)
        traj = np.transpose(np.stack(traj, axis=0), [1, 0, 2])
        return traj.squeeze()

    def _refresh_analysis_dictionary(self, sample):
        sample = tf.convert_to_tensor(sample, dtype=self.training_input_dtype)
        _ = self.dic(sample)

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

            analysis_sample = tf.convert_to_tensor(sample, dtype=self.analysis_dtype)
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
        sample = data_x[:1]
        analysis_dic = self._refresh_analysis_dictionary(sample)
        d = analysis_dic(tf.convert_to_tensor(sample, dtype=self.analysis_dtype)).shape[-1]

        # initialize accumulators
        XtX_sum = tf.zeros((d, d), dtype=self.gram_dtype)
        XtY_sum = tf.zeros((d, d), dtype=self.gram_dtype)
        n_total = 0

        total = data_x.shape[0]

        print(self.batch_size, total)
        # iterate in batches with tqdm progress bar
        for start in tqdm(range(0, total, self.batch_size), desc="compute_K batches", unit="batch"):
            end = min(start + self.batch_size, total)
            Xb = data_x[start:end]
            Yb = data_y[start:end]

            psi_x_b = analysis_dic(tf.convert_to_tensor(Xb, dtype=self.analysis_dtype))
            psi_y_b = analysis_dic(tf.convert_to_tensor(Yb, dtype=self.analysis_dtype))
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
        if data.shape[0] == 0:
            return np.empty((0, 0), dtype=self.spectral_dtype)

        analysis_dic = self._refresh_analysis_dictionary(data[:1])
        psi = analysis_dic(tf.convert_to_tensor(data, dtype=self.analysis_dtype)).numpy()
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

        for start in range(0, data_x.shape[0], batch_size):
            end = min(start + batch_size, data_x.shape[0])
            x_batch = tf.convert_to_tensor(data_x[start:end], dtype=self.training_input_dtype)
            y_batch = tf.convert_to_tensor(data_y[start:end], dtype=self.training_input_dtype)
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
        inputs_x = Input((self.target_dim,), dtype=self.training_input_dtype)
        inputs_y = Input((self.target_dim,), dtype=self.training_input_dtype)

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
    def _make_ds(
        self,
        x_array,
        y_array,
        shuffle=False,
        shuffle_buffer_size=None,
        shuffle_seed=None,
        reshuffle_each_iteration=True
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
                reshuffle_each_iteration=reshuffle_each_iteration
            )
        ds = ds.batch(self.batch_size)

        def map_to_zero(index_batch):
            def fetch_numpy_batches(index_values):
                return x_array[index_values], y_array[index_values]

            x_batch, y_batch = tf.numpy_function(
                fetch_numpy_batches,
                [index_batch],
                [self.training_input_dtype, self.training_input_dtype]
            )
            x_batch.set_shape([None, x_dim])
            y_batch.set_shape([None, y_dim])
            batch_size = tf.shape(index_batch)[0]
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
            train_shuffle=False,
            shuffle_buffer_size=None,
            shuffle_seed=None,
            reshuffle_each_iteration=True,
            outer_lr_patience=None,
            outer_lr_cooldown=0,
            outer_min_lr=0.0,
            checkpoint_path=None,
            best_checkpoint_path=None,
            save_best_only=True,
            resume=False,
            resume_mode=None,
            final_state_path=None,
            best_state_path=None,
            run_metadata=None
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
        if final_state_path is None and checkpoint_path:
            final_state_path = os.path.join(os.path.dirname(checkpoint_path), 'training_state.json')
        if best_state_path is None and best_checkpoint_path:
            best_state_path = os.path.join(os.path.dirname(best_checkpoint_path), 'training_state.json')
        final_state_path = self._normalize_writable_path(
            final_state_path,
            self.default_checkpoint_root
        )
        best_state_path = self._normalize_writable_path(
            best_state_path,
            self.default_checkpoint_root
        )
        self.final_state_path = final_state_path
        self.best_state_path = best_state_path
        run_metadata = dict(run_metadata or {})
        if resume_mode is None:
            resume_mode = 'final' if resume else 'fresh'
        resume_mode = str(resume_mode).strip().lower()
        if resume_mode not in {'fresh', 'final', 'best'}:
            raise ValueError("resume_mode must be one of: 'fresh', 'final', 'best'.")
        
        self.data_x_train, self.data_y_train = self.separate_data(self.data_train)
        print(f"build: separate_data done, x_train shape={self.data_x_train.shape}, y_train shape={self.data_y_train.shape}")
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

        self.ds_train = self._make_ds(
            self.data_x_train,
            self.data_y_train,
            shuffle=train_shuffle,
            shuffle_buffer_size=shuffle_buffer_size,
            shuffle_seed=shuffle_seed,
            reshuffle_each_iteration=reshuffle_each_iteration
        )
        self.ds_valid = self._make_ds(self.data_valid[0], self.data_valid[1])

        # self.zeros_data_y_train = tf.zeros_like(self.dic_func(self.data_y_train))
        # self.zeros_data_y_valid = tf.zeros_like(self.dic_func(self.data_valid[1]))
        print(f"creating datasets")

        self.model = self.build_model()
        print("build: build_model done")

        opt = Adam(lr)
        self.model.compile(optimizer=opt, loss=self._packed_residual_loss)
        if hasattr(opt, 'build'):
            opt.build(self.model.trainable_variables)
        self._sync_training_spectral_state()

        K_var = tf.Variable(self.K, trainable=False, name='K')
        reg_var = tf.Variable(self.reg, trainable=False, name='reg')

        checkpoint = None
        checkpoint_dir = None
        best_checkpoint_dir = None
        if checkpoint_path or best_checkpoint_path:
            if checkpoint_path:
                checkpoint_dir = os.path.dirname(checkpoint_path) or '.'
            if best_checkpoint_path:
                best_checkpoint_dir = os.path.dirname(best_checkpoint_path) or '.'
            checkpoint = tf.train.Checkpoint(
                model=self.model,
                optimizer=opt,
                K=K_var,
                eigenvectors=self.eigenvectors_var,
                eigenvalues=self.eigenvalues_var,
                reg=reg_var
            )
        resume_checkpoint_prefix = None
        resume_state_path = None
        if resume_mode == 'final':
            resume_checkpoint_prefix = checkpoint_path
            resume_state_path = final_state_path
        elif resume_mode == 'best':
            resume_checkpoint_prefix = best_checkpoint_path
            resume_state_path = best_state_path

        latest_checkpoint = None
        last_restored_checkpoint = None
        latest_final_saved_checkpoint = None
        latest_best_saved_checkpoint = None
        loaded_training_state = None
        if resume_mode == 'fresh':
            print("build: resume_mode='fresh', starting without checkpoint restore")
        elif resume_checkpoint_prefix:
            resume_checkpoint_dir = os.path.dirname(resume_checkpoint_prefix) or '.'
            latest_checkpoint = tf.train.latest_checkpoint(resume_checkpoint_dir)
            if latest_checkpoint:
                checkpoint.restore(latest_checkpoint)
                last_restored_checkpoint = latest_checkpoint
                self.K = K_var.numpy()
                self.eigenvectors = self.eigenvectors_var.numpy().astype(np.complex128, copy=False)
                self.eigenvalues = self.eigenvalues_var.numpy().astype(np.complex128, copy=False)
                self.reg = float(reg_var.numpy())
                self._sync_training_spectral_state()
                print(f"build: restored {resume_mode} checkpoint from {latest_checkpoint}")
            else:
                print(f"build: no {resume_mode} checkpoint found, starting fresh")
                resume_mode = 'fresh'

            if resume_mode != 'fresh' and resume_state_path:
                loaded_training_state = _load_json_file(resume_state_path)
                if loaded_training_state is None:
                    print(
                        f"build: no {resume_mode} trainer state found at {resume_state_path}; "
                        "continuing with checkpoint-only restore"
                    )
                else:
                    print(f"build: loaded {resume_mode} trainer state from {resume_state_path}")
        else:
            print(f"build: resume_mode='{resume_mode}' requested but no checkpoint prefix is configured; starting fresh")
            resume_mode = 'fresh'
        print(f"build: optimizer lr={opt.lr.numpy()}")

        losses = []
        val_losses = []
        learning_rate_changes = []
        stop_flag = 0
        best_val_metric = np.inf
        best_outer_epoch = None
        best_checkpoint_saved = False
        epochs_since_best = 0
        lr_cooldown_counter = 0

        if outer_lr_patience is not None:
            outer_lr_patience = int(outer_lr_patience)
            if outer_lr_patience <= 0:
                outer_lr_patience = None
        outer_lr_cooldown = max(0, int(outer_lr_cooldown))
        outer_min_lr = float(outer_min_lr)

        def _sync_checkpoint_state():
            K_var.assign(tf.cast(self.K, K_var.dtype))
            self.eigenvectors_var.assign(tf.cast(self.eigenvectors, self.training_complex_dtype))
            self.eigenvalues_var.assign(tf.cast(self.eigenvalues, self.training_complex_dtype))
            reg_var.assign(tf.cast(self.reg, reg_var.dtype))

        def _save_training_state(state_path, state_kind):
            if not state_path:
                return

            payload = {
                'schema_version': 1,
                'state_kind': state_kind,
                'resume_mode': resume_mode,
                'loss_history': list(losses),
                'val_loss_history': list(val_losses),
                'outer_history': self.outer_history,
                'best_val_metric': self.best_val_metric,
                'best_outer_epoch': self.best_outer_epoch,
                'best_outer_index': self.best_outer_index,
                'best_train_metric': self.best_train_metric,
                'best_eigvec_cond': self.best_eigvec_cond,
                'best_lr': self.best_lr,
                'best_reg': self.best_reg,
                'best_outer_summary': self.best_outer_summary,
                'best_checkpoint_saved': bool(self.best_checkpoint_saved),
                'epochs_since_best': int(epochs_since_best),
                'lr_cooldown_counter': int(lr_cooldown_counter),
                'optimizer_lr': float(opt.lr.numpy()),
                'current_reg': float(self.reg),
                'current_eigvec_cond': float(self.eigvec_cond) if self.eigvec_cond is not None else None,
                'completed_outer_epochs': int(len(self.outer_history)),
                'completed_rounds': int(run_metadata.get('completed_rounds', 0)),
                'round_number': int(run_metadata.get('round_number', 0)),
                'run_metadata': run_metadata,
                'loaded_checkpoint': latest_checkpoint,
                'restored_checkpoint': last_restored_checkpoint,
                'checkpoint_prefixes': {
                    'final': checkpoint_path,
                    'best': best_checkpoint_path,
                },
                'latest_saved_checkpoints': {
                    'final': latest_final_saved_checkpoint,
                    'best': latest_best_saved_checkpoint,
                },
                'state_paths': {
                    'final': final_state_path,
                    'best': best_state_path,
                },
            }
            _save_json_file(state_path, payload)

        def _safe_float(value, default=np.nan):
            try:
                if value is None:
                    return float(default)
                return float(value)
            except Exception:
                return float(default)

        if loaded_training_state is not None:
            losses = [float(x) for x in loaded_training_state.get('loss_history', [])]
            val_losses = [float(x) for x in loaded_training_state.get('val_loss_history', [])]
            self.outer_history = list(loaded_training_state.get('outer_history', []))
            self.best_val_metric = _safe_float(loaded_training_state.get('best_val_metric', np.nan))
            self.best_outer_epoch = loaded_training_state.get('best_outer_epoch', None)
            self.best_outer_index = loaded_training_state.get('best_outer_index', None)
            self.best_train_metric = _safe_float(loaded_training_state.get('best_train_metric', np.nan))
            self.best_eigvec_cond = _safe_float(loaded_training_state.get('best_eigvec_cond', np.nan))
            self.best_lr = _safe_float(loaded_training_state.get('best_lr', np.nan))
            self.best_reg = _safe_float(loaded_training_state.get('best_reg', np.nan))
            self.best_outer_summary = dict(loaded_training_state.get('best_outer_summary', {}) or {})
            self.best_checkpoint_saved = bool(
                loaded_training_state.get('best_checkpoint_saved', False)
            )
            best_val_metric = (
                self.best_val_metric if np.isfinite(self.best_val_metric) else np.inf
            )
            best_outer_epoch = self.best_outer_epoch
            best_checkpoint_saved = self.best_checkpoint_saved
            epochs_since_best = int(loaded_training_state.get('epochs_since_best', 0))
            lr_cooldown_counter = int(loaded_training_state.get('lr_cooldown_counter', 0))
            restored_lr = loaded_training_state.get('optimizer_lr', None)
            if restored_lr is not None and np.isfinite(restored_lr):
                opt.lr.assign(float(restored_lr))
            print(
                "build: restored trainer state "
                f"(outer_history={len(self.outer_history)}, loss_len={len(losses)}, "
                f"val_loss_len={len(val_losses)})"
            )

        outer_epoch_offset = len(self.outer_history)
        target_outer_epoch = outer_epoch_offset + int(epochs)

        for i in range(epochs):
            outer_epoch_start = time.time()
            current_outer_epoch = outer_epoch_offset + i + 1
            print(f"build: Outer Epoch {current_outer_epoch}/{target_outer_epoch}")
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
                callbacks=[EarlyStopping(monitor='val_loss', min_delta=1e-6, patience=2, verbose=1, mode='min', restore_best_weights=True)]
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
                'outer_epoch': current_outer_epoch,
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
                best_outer_epoch = current_outer_epoch
                self.best_val_metric = float(val_metric)
                self.best_outer_epoch = current_outer_epoch
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
                epochs_since_best = 0
            else:
                epochs_since_best += 1

            if i % log_interval == 0:
                losses.extend(self.history.history['loss'])
                val_losses.extend(self.history.history['val_loss'])

            should_save_best_checkpoint = (
                checkpoint is not None
                and best_checkpoint_path is not None
                and (current_is_better or not save_best_only)
            )
            if should_save_best_checkpoint:
                if best_checkpoint_dir and not os.path.exists(best_checkpoint_dir):
                    os.makedirs(best_checkpoint_dir, exist_ok=True)
                _sync_checkpoint_state()
                latest_best_saved_checkpoint = checkpoint.save(file_prefix=best_checkpoint_path)
                best_checkpoint_saved = True
                self.best_checkpoint_saved = True
                if current_is_better:
                    print("build: best checkpoint saved")
                _save_training_state(best_state_path, 'best')

            if lr_cooldown_counter > 0:
                lr_cooldown_counter -= 1
            elif outer_lr_patience is not None and epochs_since_best >= outer_lr_patience:
                old_lr = float(opt.lr.numpy())
                new_lr = max(old_lr * lr_decay_factor, outer_min_lr)
                if new_lr < old_lr:
                    print(
                        "build: no outer best improvement for "
                        f"{epochs_since_best} epoch(s), decaying lr "
                        f"from {old_lr:.6e} to {new_lr:.6e}"
                    )
                    opt.lr.assign(new_lr)
                else:
                    print(
                        "build: outer lr plateau reached but lr is already at "
                        f"the configured minimum ({outer_min_lr:.6e})"
                    )
                epochs_since_best = 0
                lr_cooldown_counter = outer_lr_cooldown

            learning_rate_changes.append(opt.lr.numpy())
            _save_training_state(final_state_path, 'final')
            print(f"build: Outer Epoch {current_outer_epoch} finished in {time.time() - outer_epoch_start:.2f}s")

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
            latest_final_saved_checkpoint = checkpoint.save(file_prefix=checkpoint_path)
            print("build: checkpoint saved")
        _save_training_state(final_state_path, 'final')

        return losses, val_losses, stop_flag, learning_rate_changes
