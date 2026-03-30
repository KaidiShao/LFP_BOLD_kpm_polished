print("edmd_utils module loaded")
import tensorflow as tf
import numpy as np
import sys
import os
import scipy.io
import h5py
import importlib

# Function to set device configuration
def set_device(device_choice):
    if device_choice.lower() == 'gpu':
        gpus = tf.config.experimental.list_physical_devices('GPU')
        if gpus:
            try:
                for gpu in gpus:
                    tf.config.experimental.set_memory_growth(gpu, True)
                print("Using GPU")
            except RuntimeError as e:
                print(e)
        else:
            print("GPU not available, using CPU")
    elif device_choice.lower() == 'cpu':
        cpus = tf.config.experimental.list_physical_devices('CPU')
        tf.config.set_visible_devices([], 'GPU')
        tf.config.set_visible_devices(cpus, 'CPU')
        print("Using CPU")


def import_solver(solver_name, search_dirs=None):
    # search_dirs can be overridden from a notebook/script if needed
    if search_dirs is None:
        search_dirs = [
            '/root/projects/solver',
            '/root/projects'
        ]

    for d in search_dirs:
        if os.path.isdir(d) and d not in sys.path:
            sys.path.append(d)

    module_candidates = {
        'edmd': ['solver_edmd'],
        'resdmd': ['solver_resdmd'],
        'resdmd_batch': ['solver_resdmd_batch', 'solver_resdmd_batch2'],
        'resdmd2': ['solver_resdmd2'],
        'resdmd_torch': ['solver_resdmd_torch_gpu7']
    }

    key = solver_name.lower()
    if key not in module_candidates:
        raise ValueError(f"Unknown solver name: {solver_name}")

    last_error = None
    for mod_name in module_candidates[key]:
        try:
            if mod_name in sys.modules:
                mod = importlib.reload(sys.modules[mod_name])
            else:
                mod = importlib.import_module(mod_name)
            return mod.KoopmanSolver
        except ModuleNotFoundError as e:
            last_error = e

    raise ModuleNotFoundError(
        f"Could not import solver '{solver_name}'. Tried {module_candidates[key]} "
        f"from search_dirs={search_dirs}"
    ) from last_error


def load_data(filename, file_path, file_type, field_name):
    full_path = os.path.join(file_path, filename)
    if file_type.lower() == '.mat':
        mat = scipy.io.loadmat(full_path)
        loaded_data = mat[field_name]
    elif file_type.lower() == '.h5':
        with h5py.File(full_path, 'r') as file:
            loaded_data = np.array(file['/' + field_name])
    else:
        raise ValueError(f"Unsupported file_type: {file_type}")

    print('size of loaded data: ' + str(loaded_data.shape))
    return loaded_data


def load_data_from_XY(filename, file_path, file_type, field_name_x, field_name_y):
    full_path = os.path.join(file_path, filename)
    if file_type.lower() == '.mat':
        mat = scipy.io.loadmat(full_path)
        loaded_data_x = mat[field_name_x]
        loaded_data_y = mat[field_name_y]
    elif file_type.lower() == '.h5':
        with h5py.File(full_path, 'r') as file:
            loaded_data_x = np.array(file['/' + field_name_x])
            loaded_data_y = np.array(file['/' + field_name_y])
    else:
        raise ValueError(f"Unsupported file_type: {file_type}")

    print('size of loaded data_x: ' + str(loaded_data_x.shape))
    print('size of loaded data_y: ' + str(loaded_data_y.shape))
    return loaded_data_x, loaded_data_y


def remove_checkpoint(checkpoint_filepath):
    if os.path.exists(checkpoint_filepath):
        import shutil
        shutil.rmtree(checkpoint_filepath)
        print(f"Removed checkpoint directory: {checkpoint_filepath}")
    else:
        print(f"No checkpoint directory found at: {checkpoint_filepath}")


def transfer_data_format(data, train_ratio=0.7):
    X = data[0:-1, :]
    Y = data[1:, :]

    len_all = X.shape[0]

    np.random.seed(100)
    random_indices = np.random.permutation(len_all)

    train_indices = random_indices[:int(train_ratio * len_all)]
    test_indices = random_indices[int(train_ratio * len_all):]

    data_x_train = X[train_indices, :]
    data_x_valid = X[test_indices, :]

    data_y_train = Y[train_indices, :]
    data_y_valid = Y[test_indices, :]

    data_train = [data_x_train, data_y_train]
    data_valid = [data_x_valid, data_y_valid]

    return data_train, data_valid


def dim_reduction(data, n_modes):
    U, S, V = np.linalg.svd(data, full_matrices=False)

    U = U[:, 0:n_modes]
    S = np.diag(S[0:n_modes])
    V = V[0:n_modes, :]

    data_reduced = np.dot(U, S)

    return data_reduced, U, S, V
