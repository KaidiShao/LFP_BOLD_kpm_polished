import os
# from autograd import jacobian, hessian
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
from tqdm import tqdm


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.keras.backend.set_floatx('float32')

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
            spectral_dtype='float64'):
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
        self.training_policy = training_policy or (
            'mixed_float16' if tf.config.list_physical_devices('GPU') else 'float32'
        )
        self.analysis_dtype = tf.as_dtype(analysis_dtype)
        self.gram_dtype = tf.as_dtype(gram_dtype)
        self.spectral_dtype = np.dtype(spectral_dtype)
        self.training_input_dtype = tf.float32
        self.psi_x = None
        self.psi_y = None
        self.analysis_dic = None
        self.output_dim = None
        self.eigenvectors_var = None
        mixed_precision.set_global_policy(self.training_policy)

    def separate_data(self, data):
        data_x = data[0]
        data_y = data[1]
        return data_x, data_y

    def build(self, data_train):
        # Separate data
        self.data_train = data_train
        self.data_x_train, self.data_y_train = self.separate_data(
            self.data_train)

        # Compute final information
        self.compute_final_info(reg_final=0.0)

    def compute_final_info(self, reg_final):
        # Compute K
        self.K = self.compute_K(self.dic_func,
                                self.data_x_train,
                                self.data_y_train,
                                reg=reg_final)
        self.eig_decomp(self.K)
        self.compute_mode()

    def eig_decomp(self, K):
        """ eigen-decomp of K """
        K_np = np.asarray(K, dtype=self.spectral_dtype)
        self.eigenvalues, self.eigenvectors = np.linalg.eig(K_np)
        idx = self.eigenvalues.real.argsort()[::-1]
        self.eigenvalues = self.eigenvalues[idx]
        self.eigenvectors = self.eigenvectors[:, idx]
        self.eigenvectors_inv = np.linalg.inv(self.eigenvectors)

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
        self.modes = np.matmul(self.eigenvectors_inv, self.B).T
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
        sample = data_x[:1]
        analysis_dic = self._refresh_analysis_dictionary(sample)
        d = analysis_dic(tf.convert_to_tensor(sample, dtype=self.analysis_dtype)).shape[-1]

        # initialize accumulators
        XtX = tf.zeros((d, d), dtype=self.gram_dtype)
        XtY = tf.zeros((d, d), dtype=self.gram_dtype)

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

            # accumulate Gram and cross-covariance
            XtX += tf.matmul(psi_x_b, psi_x_b, transpose_a=True)
            XtY += tf.matmul(psi_x_b, psi_y_b, transpose_a=True)

        # form and solve normal equations
        idmat = tf.eye(d, dtype=self.gram_dtype)
        inv_term = tf.linalg.pinv(XtX + reg * idmat)
        K_reg = tf.matmul(inv_term, XtY)

        return K_reg




    def get_Psi_X(self):
        return self.Psi_X

    def get_Psi_Y(self):
        return self.Psi_Y


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

        Layer_K = Dense(units=self.psi_y.shape[-1],
                        use_bias=False,
                        name='Layer_K',
                        trainable=False,
                        dtype='float32')
        psi_next = Layer_K(tf.cast(self.psi_x, tf.float32))

        if self.eigenvectors_var is None:
            raise ValueError("self.eigenvectors_var must be initialized before build_model().")

        psi_diff = tf.cast(psi_next - tf.cast(self.psi_y, tf.float32), tf.complex64)
        weighted_diff = tf.matmul(psi_diff, tf.cast(self.eigenvectors_var, tf.complex64))
        outputs = tf.concat(
            [tf.math.real(weighted_diff), tf.math.imag(weighted_diff)],
            axis=-1
        )
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
        x_tensor = tf.convert_to_tensor(x_array, dtype=self.training_input_dtype)
        y_tensor = tf.convert_to_tensor(y_array, dtype=self.training_input_dtype)
        ds = tf.data.Dataset.from_tensor_slices((x_tensor, y_tensor))
        ds = ds.batch(self.batch_size)

        def map_to_zero(x_batch, y_batch):
            batch_size = tf.shape(x_batch)[0]
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
            end_condition=1e-5
        ):
        print("build: start")
        self.data_train = data_train
        self.batch_size = batch_size
        
        self.data_x_train, self.data_y_train = self.separate_data(self.data_train)
        print(f"build: separate_data done, x_train shape={self.data_x_train.shape}, y_train shape={self.data_y_train.shape}")
        self.data_valid = data_valid

        print("build: compute_final_info start")
        self.compute_final_info(reg_final=0.01)
        print(f"build: compute_final_info done, K shape={self.K.shape}")
        self.output_dim = int(self.K.shape[0]) * 2
        self.eigenvectors_var = tf.Variable(self.eigenvectors, trainable=False, name='eigenvectors')

        self.ds_train = self._make_ds(self.data_x_train, self.data_y_train)
        self.ds_valid = self._make_ds(self.data_valid[0], self.data_valid[1])

        # self.zeros_data_y_train = tf.zeros_like(self.dic_func(self.data_y_train))
        # self.zeros_data_y_valid = tf.zeros_like(self.dic_func(self.data_valid[1]))
        print(f"creating datasets")

        self.model = self.build_model()
        print("build: build_model done")

        opt = Adam(lr)
        self.model.compile(optimizer=opt, loss='mse')

        checkpoint_path = './checkpoints_test/'
        checkpoint_dir = os.path.dirname(checkpoint_path)
        if not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)

        K_var = tf.Variable(self.K, name='K')
        eigenvalues_var = tf.Variable(self.eigenvalues, name='eigenvalues')
        reg_var = tf.Variable(self.reg, name='reg')

        dummy_gradients = [tf.zeros_like(var) for var in self.model.trainable_variables]
        opt.apply_gradients(zip(dummy_gradients, self.model.trainable_variables))

        checkpoint = tf.train.Checkpoint(
            model=self.model,
            optimizer=opt,
            K=K_var,
            eigenvectors=self.eigenvectors_var,
            eigenvalues=eigenvalues_var,
            reg=reg_var
        )
        latest_checkpoint = tf.train.latest_checkpoint(checkpoint_dir)
        if latest_checkpoint:
            checkpoint.restore(latest_checkpoint)
            self.K = K_var.numpy()
            self.eigenvectors = self.eigenvectors_var.numpy()
            self.eigenvalues = eigenvalues_var.numpy()
            self.reg = reg_var.numpy()
            print(f"build: restored checkpoint from {latest_checkpoint}")
        else:
            print("build: no checkpoint found, starting fresh")
        print(f"build: optimizer lr={opt.lr.numpy()}")

        losses = []
        val_losses = []
        learning_rate_changes = []
        stop_flag = 0

        for i in range(epochs):
            print(f"build: Outer Epoch {i+1}/{epochs}")
            print("build: compute_K start")
            self.K = self.compute_K(self.dic_func, self.data_x_train, self.data_y_train, self.reg)
            self.eig_decomp(self.K)
            self.eigenvectors_var.assign(self.eigenvectors)
            print("build: compute_K done")

            layer_k_weights = self.model.get_layer('Layer_K').weights[0]
            layer_k_weights.assign(tf.cast(self.K, layer_k_weights.dtype))
            print("build: assigned K to Layer_K")

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

            learning_rate_changes.append(opt.lr.numpy())

            if i % log_interval == 0:
                losses.extend(self.history.history['loss'])
                val_losses.extend(self.history.history['val_loss'])
                if len(losses) > 2 and losses[-1] > losses[-2]:
                    print("build: loss increased, decaying lr")
                    opt.lr.assign(opt.lr * lr_decay_factor)
                if len(val_losses) > 200:
                    avg_loss_change = np.mean([abs(val_losses[j] - val_losses[-1]) for j in range(-200, 0)])
                    if avg_loss_change < end_condition:
                        print("build: loss stabilized, stopping training")
                        stop_flag = 1
                        break

        print("build: final compute_final_info start")
        self.compute_final_info(reg_final=0.01)
        print("build: final compute_final_info done")

        layer_k_weights = self.model.get_layer('Layer_K').weights[0]
        layer_k_weights.assign(tf.cast(self.K, layer_k_weights.dtype))
        K_var.assign(self.K)
        self.eigenvectors_var.assign(self.eigenvectors)
        eigenvalues_var.assign(self.eigenvalues)
        reg_var.assign(self.reg)

        opt.apply_gradients(zip(dummy_gradients, self.model.trainable_variables))
        checkpoint.save(file_prefix=checkpoint_path)
        print("build: checkpoint saved")

        return losses, val_losses, stop_flag, learning_rate_changes
