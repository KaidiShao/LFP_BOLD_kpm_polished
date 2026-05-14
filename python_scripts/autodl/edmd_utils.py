print("edmd_utils module loaded")
import tensorflow as tf
import numpy as np
import sys
import os
import scipy.io
import h5py
import importlib

# Function to set device configuration
def set_device(device_choice, require_gpu=False):
    if device_choice.lower() == 'gpu':
        gpus = tf.config.experimental.list_physical_devices('GPU')
        if gpus:
            try:
                for gpu in gpus:
                    tf.config.experimental.set_memory_growth(gpu, True)
                print(f"Using GPU: {gpus}")
                return 'gpu'
            except RuntimeError as e:
                print(e)
        else:
            if require_gpu:
                raise RuntimeError("GPU was requested with require_gpu=True, but TensorFlow found no GPU devices.")
            print("GPU not available, using CPU")
            return 'cpu'
    elif device_choice.lower() == 'cpu':
        cpus = tf.config.experimental.list_physical_devices('CPU')
        tf.config.set_visible_devices([], 'GPU')
        tf.config.set_visible_devices(cpus, 'CPU')
        print("Using CPU")
        return 'cpu'
    return device_choice.lower()


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
        'resdmd_batch': ['solver_resdmd_batch3', 'solver_resdmd_batch2', 'solver_resdmd_batch'],
        'resdmd_batch4': ['solver_resdmd_batch4'],
        'resdmd_batch_static': ['solver_resdmd_batch3_static_target'],
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


def _hdf5_dataset_key(file, field_name):
    key = field_name if field_name.startswith('/') else '/' + field_name
    if key not in file:
        available = ', '.join(sorted(file.keys()))
        raise KeyError(f"Dataset {key} not found in HDF5 file. Available top-level keys: {available}")
    return key


def _load_hdf5_dataset(full_path, field_name):
    with h5py.File(full_path, 'r') as file:
        key = _hdf5_dataset_key(file, field_name)
        return np.array(file[key])


def load_data(filename, file_path, file_type, field_name):
    full_path = os.path.join(file_path, filename)
    if file_type.lower() == '.mat':
        if h5py.is_hdf5(full_path):
            loaded_data = _load_hdf5_dataset(full_path, field_name)
        else:
            mat = scipy.io.loadmat(full_path)
            loaded_data = mat[field_name]
    elif file_type.lower() == '.h5':
        loaded_data = _load_hdf5_dataset(full_path, field_name)
    else:
        raise ValueError(f"Unsupported file_type: {file_type}")

    print('size of loaded data: ' + str(loaded_data.shape))
    return loaded_data


def load_data_from_XY(filename, file_path, file_type, field_name_x, field_name_y):
    full_path = os.path.join(file_path, filename)
    if file_type.lower() == '.mat':
        if h5py.is_hdf5(full_path):
            loaded_data_x = _load_hdf5_dataset(full_path, field_name_x)
            loaded_data_y = _load_hdf5_dataset(full_path, field_name_y)
        else:
            mat = scipy.io.loadmat(full_path)
            loaded_data_x = mat[field_name_x]
            loaded_data_y = mat[field_name_y]
    elif file_type.lower() == '.h5':
        loaded_data_x = _load_hdf5_dataset(full_path, field_name_x)
        loaded_data_y = _load_hdf5_dataset(full_path, field_name_y)
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


class IndexedRowsView:
    """Array-like view that materializes rows from a base matrix on demand."""

    def __init__(self, data, row_indices, name=None):
        self.data = data
        self.row_indices = np.asarray(row_indices, dtype=np.int64).reshape(-1)
        self.name = str(name or "indexed_rows")
        self.shape = (int(self.row_indices.size), int(data.shape[1]))
        self.dtype = getattr(data, "dtype", None)

    def __len__(self):
        return self.shape[0]

    def __getitem__(self, key):
        if isinstance(key, tuple):
            row_key, col_key = key
        else:
            row_key, col_key = key, slice(None)
        resolved_rows = self.row_indices[row_key]
        return self.data[resolved_rows, col_key]

    def get_rows(self, local_indices):
        local_indices = np.asarray(local_indices, dtype=np.int64).reshape(-1)
        return self.data[self.row_indices[local_indices], :]


def transfer_data_format(data, train_ratio=0.7, seed=100):
    len_all = int(data.shape[0] - 1)
    valid_idx_x = np.arange(len_all, dtype=np.int64)
    valid_idx_y = valid_idx_x + 1

    rng = np.random.default_rng(seed)
    random_indices = rng.permutation(len_all)

    train_indices = random_indices[:int(train_ratio * len_all)]
    test_indices = random_indices[int(train_ratio * len_all):]

    data_x_train = IndexedRowsView(data, valid_idx_x[train_indices], name="train_x")
    data_x_valid = IndexedRowsView(data, valid_idx_x[test_indices], name="valid_x")

    data_y_train = IndexedRowsView(data, valid_idx_y[train_indices], name="train_y")
    data_y_valid = IndexedRowsView(data, valid_idx_y[test_indices], name="valid_y")

    data_train = [data_x_train, data_y_train]
    data_valid = [data_x_valid, data_y_valid]

    return data_train, data_valid


def make_session_aware_snapshot_indices(session_start_idx, session_end_idx, lag=1, matlab_indexing=True):
    starts = np.asarray(session_start_idx).reshape(-1).astype(np.int64)
    ends = np.asarray(session_end_idx).reshape(-1).astype(np.int64)

    if starts.size != ends.size:
        raise ValueError("session_start_idx and session_end_idx must have the same length")
    if lag < 1:
        raise ValueError("lag must be a positive integer")

    if matlab_indexing:
        starts = starts - 1
        ends = ends - 1

    valid_x = []
    valid_y = []
    session_idx = []

    for i, (start, end) in enumerate(zip(starts, ends)):
        if end - start < lag:
            continue
        ix = np.arange(start, end - lag + 1, dtype=np.int64)
        iy = ix + lag
        valid_x.append(ix)
        valid_y.append(iy)
        session_idx.append(np.full(ix.shape, i, dtype=np.int64))

    if not valid_x:
        raise ValueError("No valid session-aware snapshot pairs were generated")

    return (
        np.concatenate(valid_x),
        np.concatenate(valid_y),
        np.concatenate(session_idx),
    )


def transfer_data_format_session_aware(
    data,
    session_start_idx,
    session_end_idx,
    train_ratio=0.7,
    lag=1,
    seed=100,
    matlab_indexing=True,
):
    valid_idx_x, valid_idx_y, session_idx = make_session_aware_snapshot_indices(
        session_start_idx,
        session_end_idx,
        lag=lag,
        matlab_indexing=matlab_indexing,
    )

    len_all = int(valid_idx_x.shape[0])

    rng = np.random.default_rng(seed)
    random_indices = rng.permutation(len_all)

    split_idx = int(train_ratio * len_all)
    train_indices = random_indices[:split_idx]
    test_indices = random_indices[split_idx:]

    data_train = [
        IndexedRowsView(data, valid_idx_x[train_indices], name="train_x"),
        IndexedRowsView(data, valid_idx_y[train_indices], name="train_y"),
    ]
    data_valid = [
        IndexedRowsView(data, valid_idx_x[test_indices], name="valid_x"),
        IndexedRowsView(data, valid_idx_y[test_indices], name="valid_y"),
    ]

    snapshot_metadata = {
        "valid_idx_x": valid_idx_x,
        "valid_idx_y": valid_idx_y,
        "session_idx": session_idx,
        "train_pair_idx": train_indices,
        "valid_pair_idx": test_indices,
        "lag": int(lag),
    }

    return data_train, data_valid, snapshot_metadata


def dim_reduction(data, n_modes):
    U, S, V = np.linalg.svd(data, full_matrices=False)

    U = U[:, 0:n_modes]
    S = np.diag(S[0:n_modes])
    V = V[0:n_modes, :]

    data_reduced = np.dot(U, S)

    return data_reduced, U, S, V
