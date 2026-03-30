print("edmd_utils module loaded")
import tensorflow as tf
import numpy as np
import sys
import os
import scipy.io
import h5py

# Function to set device configuration
def set_device(device_choice):
    if device_choice.lower() == 'gpu':
        gpus = tf.config.experimental.list_physical_devices('GPU')
        if gpus:
            try:
                # Set memory growth if using GPU
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
    

# import selected edmd solver
# def import_solver(solver_name):
#     # Ensure the path to the solver modules is in sys.path
#     sys.path.append('/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/code/algorithms/new_kpm_code_yuanchao/')
#     sys.path.append('/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/code/algorithms/new_kpm_code_yuanchao_gpu20240702/')
   
    
#     if solver_name.lower() == 'edmd':
#         import solver_edmd
#         from solver_edmd import KoopmanSolver
#     elif solver_name.lower() == 'resdmd':
#         import solver_resdmd
#         from solver_resdmd import KoopmanSolver
#     elif solver_name.lower() == 'resdmd_batch':
#         import solver_resdmd_batch
#         from solver_resdmd_batch import KoopmanSolver
#     elif solver_name.lower() == 'resdmd2':
#         import solver_resdmd2
#         from solver_resdmd2 import KoopmanSolver
#     elif solver_name.lower() == 'resdmd_torch':
#         import solver_resdmd_torch_gpu7
#         from solver_resdmd_torch_gpu7 import KoopmanSolver
#     else:
#         raise ValueError(f"Unknown solver name: {solver_name}")
    
#     return KoopmanSolver
import sys
import importlib

def import_solver(solver_name):
    # 1) make sure your code directories are on the path
    base_dirs = [
        '/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/code/algorithms/new_kpm_code_yuanchao/',
        '/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/code/algorithms/new_kpm_code_yuanchao_gpu20240702/'
    ]
    for d in base_dirs:
        if d not in sys.path:
            sys.path.append(d)
    
    # 2) list of all solver module names
    module_names = {
        'edmd':            'solver_edmd',
        'resdmd':         'solver_resdmd',
        'resdmd_batch':   'solver_resdmd_batch2',
        'resdmd2':         'solver_resdmd2',
        'resdmd_torch':   'solver_resdmd_torch_gpu7'
    }
    
    key = solver_name.lower()
    if key not in module_names:
        raise ValueError(f"Unknown solver name: {solver_name}")
    
    mod_name = module_names[key]
    
    # 3) import (or re-import) each module so any edits on disk take effect
    #    we reload all of them to clear caches
    for m in module_names.values():
        if m in sys.modules:
            importlib.reload(sys.modules[m])
        else:
            importlib.import_module(m)
    
    # 4) now reload the one we actually need (just in case)
    mod = importlib.reload(sys.modules[mod_name])
    
    # 5) return its KoopmanSolver
    return mod.KoopmanSolver




def load_data(filename, file_path, file_type, field_name):
    # note: filename should include file extention (because sometimes h5 file also has suffix of .mat)
    if file_type.lower() == '.mat':
        mat = scipy.io.loadmat(file_path + filename)
        loaded_data = mat[field_name]
    if file_type.lower() == '.h5':
        with h5py.File(file_path +  filename, 'r') as file:
            loaded_data = np.array(file['/' + field_name])  # Transpose the data

    print('size of loaded data: ' + str(loaded_data.shape))
    return loaded_data

def load_data_from_XY(filename, file_path):
    """
    Loads the variables X and Y from a .mat file.

    Parameters:
        filename (str): The name of the .mat file (including extension).
        file_path (str): The directory path where the file is located.
    
    Returns:
        tuple: A tuple (X, Y) containing the loaded variables.
    """
    full_path = os.path.join(file_path, filename)
    mat = scipy.io.loadmat(full_path)
    X = mat['X']
    Y = mat['Y']
    print(f"Loaded X with shape {X.shape} and Y with shape {Y.shape}")
    return X, Y

# remove checkpoint directory at the beginning of all rounds
def remove_checkpoint(checkpoint_filepath):
    if os.path.exists(checkpoint_filepath):
        import shutil
        shutil.rmtree(checkpoint_filepath)
        print(f"Removed checkpoint directory: {checkpoint_filepath}")
    else:
        print(f"No checkpoint directory found at: {checkpoint_filepath}")

# Function to transfer data format (divide data into training and testing sets)
def transfer_data_format(data, train_ratio=0.7):
    # data: the data to be transferred: make sure dimensions are (n_features, n_samples)
    # train_ratio: the ratio of data to be used for training
    # return: training data and testing data

    # define snapshots
    X = data[0:-1,:]
    Y = data[1:,:]

    # define training and testing indices
    len_all = X.shape[0]

    np.random.seed(100)
    random_indices = np.random.permutation(len_all)

    train_indices = random_indices[:int(0.7 * len_all)]
    test_indices = random_indices[int(0.7 * len_all):]

    # extract training and testing data
    data_x_train = X[train_indices,:]
    data_x_valid = X[test_indices,:]

    data_y_train = Y[train_indices,:]
    data_y_valid = Y[test_indices,:]

    data_train = [data_x_train, data_y_train]
    data_valid = [data_x_valid, data_y_valid]

    return data_train, data_valid

def transfer_data_format_from_XY(X, Y, train_ratio=0.7):
    # data: the data to be transferred: make sure dimensions are (n_features, n_samples)
    # train_ratio: the ratio of data to be used for training
    # return: training data and testing data

    # define snapshots
    X = X.T
    Y = Y.T

    # define training and testing indices
    len_all = X.shape[0]

    random_indices = np.random.permutation(len_all)

    train_indices = random_indices[:int(0.7 * len_all)]
    test_indices = random_indices[int(0.7 * len_all):]

    # extract training and testing data
    data_x_train = X[train_indices,:]
    data_x_valid = X[test_indices,:]

    data_y_train = Y[train_indices,:]
    data_y_valid = Y[test_indices,:]

    data_train = [data_x_train, data_y_train]
    data_valid = [data_x_valid, data_y_valid]

    return data_train, data_valid

def transfer_data_format_sensorium(data, train_ratio=0.7):
    # data: the data to be transferred: make sure dimensions are (n_features, n_samples)
    # train_ratio: the ratio of data to be used for training
    # return: training data and testing data

    # define snapshots
    # X = data[0:-1,:]
    # Y = data[1:,:]
    # Step 1: Calculate the first dimension and divide by 300
    n_segments = data.shape[0] // 300  # Calculate 58 as n_segments 

    # Step 2: Reshape to (n_segments, 300, 7863)
    reshaped_data = data.reshape(n_segments, 300, data.shape[1])

    # Step 3: Extract indices 1:299 from the second dimension
    X = reshaped_data[:, 0:299, :].reshape(-1, data.shape[1])  # Shape will be (n_segments, 299, 7863)
    Y = reshaped_data[:, 1:300, :].reshape(-1, data.shape[1])  # Shape will be (n_segments, 299, 7863)
    print(X.shape, Y.shape)

    # define training and testing indices
    len_all = X.shape[0]

    random_indices = np.random.permutation(len_all)

    train_indices = random_indices[:int(0.7 * len_all)]
    test_indices = random_indices[int(0.7 * len_all):]

    # extract training and testing data
    data_x_train = X[train_indices,:]
    data_x_valid = X[test_indices,:]

    data_y_train = Y[train_indices,:]
    data_y_valid = Y[test_indices,:]

    data_train = [data_x_train, data_y_train]
    data_valid = [data_x_valid, data_y_valid]

    return data_train, data_valid

def dim_reduction(data, n_modes):
    # data: the data to be reduced: make sure dimensions are (n_features, n_samples)
    # n_modes: the number of modes to be retained
    # return: reduced data

    # SVD
    U, S, V = np.linalg.svd(data, full_matrices=False)

    # retain the first n_modes
    U = U[:,0:n_modes]
    S = np.diag(S[0:n_modes])
    V = V[0:n_modes,:]

    # reduced data
    data_reduced = np.dot(U, S)

    return data_reduced, U, S, V
