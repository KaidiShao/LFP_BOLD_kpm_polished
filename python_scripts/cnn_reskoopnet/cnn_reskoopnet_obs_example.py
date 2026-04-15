import numpy as np

from cnn_reskoopnet_data import build_obs_1dcnn_dataset, load_obs_info_csv
from solver_resdmd_cnn import KoopmanAuxCNN, KoopmanSolver

try:
    import h5py
except ImportError:  # pragma: no cover - depends on local environment
    h5py = None


def load_v73_mat_field(mat_path, field_name):
    if h5py is None:
        raise ImportError(
            "h5py is required to read MATLAB v7.3/HDF5 files in this example."
        )

    with h5py.File(mat_path, 'r') as f:
        return np.array(f['/' + field_name])


def load_obs_and_border_idx(mat_path, obs_field='obs', border_field='border_idx'):
    obs = load_v73_mat_field(mat_path, obs_field).T.astype(np.float32, copy=False)
    border_idx = load_v73_mat_field(mat_path, border_field).reshape(-1).astype(int, copy=False)
    return obs, border_idx


if __name__ == '__main__':
    # Current F12m01 dictionary files:
    # - abs: spec becomes [time, 151, 2] = [freq_features, hp/pl]
    # - complex_split: spec becomes [time, 151, 4] = [hp_real, hp_imag, pl_real, pl_imag]
    data_path = r'D:\DataPons_processed\f12m01\reskoopnet_dictionary\f12m01_low50_high250_g2_abs_single.mat'
    info_csv_path = r'D:\DataPons_processed\f12m01\reskoopnet_dictionary\f12m01_low50_high250_g2_abs_single_obs_info.csv'

    obs, border_idx = load_obs_and_border_idx(data_path)
    obs_info_rows = load_obs_info_csv(info_csv_path)

    data_train, data_valid, metadata = build_obs_1dcnn_dataset(
        obs,
        obs_info_rows,
        border_idx=border_idx,
        train_ratio=0.7,
        seed=100,
        output_structure='dict'
    )

    print('dataset metadata:', metadata)

    dic = KoopmanAuxCNN(
        spec_input_shape=(metadata['spec_n_freq_features'], metadata['spec_n_channels']),
        aux_dim=metadata['aux_dim'],
        conv_filters=(16, 32, 64),
        kernel_size=3,
        pool_size=2,
        dense_units=(64,),
        n_psi_train=32,
    )

    solver = KoopmanSolver(
        dic=dic,
        target_dim=metadata['aux_dim'],
        reg=1e-3,
    )

    losses, val_losses, stop_flag, learning_rate_changes = solver.build(
        data_train=data_train,
        data_valid=data_valid,
        epochs=10,
        batch_size=512,
        lr=1e-3,
        log_interval=1,
        lr_decay_factor=0.5,
        Nepoch=5,
    )

    print('final stop_flag:', stop_flag)
    print('num losses:', len(losses), 'num val losses:', len(val_losses))
    print('lr history length:', len(learning_rate_changes))
