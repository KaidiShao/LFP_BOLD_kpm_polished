import argparse
import gc
import os
import shutil
import sys
from pathlib import Path


def resolve_local_results_root():
    env_root = os.environ.get('LFP_BOLD_KPM_RESULTS_ROOT')
    if env_root:
        return Path(env_root)
    if os.name == 'nt':
        return Path('E:/autodl_results')
    return Path('/mnt/e/autodl_results')


def sanitize_tag(value):
    return str(value).replace('/', '-').replace(' ', '_')


def parse_args():
    default_data_root = Path('/mnt/d/DataPons_processed/f12m01/reskoopnet_dictionary')

    parser = argparse.ArgumentParser(
        description='Run the CNN-based ResKoopNet pipeline on the full obs dataset.'
    )
    parser.add_argument('--dataset-stem', default='f12m01')
    parser.add_argument('--observable-mode', default='abs', choices=['abs', 'complex_split'])
    parser.add_argument('--device', default='cpu', choices=['cpu', 'gpu'])
    parser.add_argument('--data-root', type=Path, default=default_data_root)
    parser.add_argument('--data-filename', default=None)
    parser.add_argument('--info-csv-filename', default=None)
    parser.add_argument('--output-parent', type=Path, default=None)
    parser.add_argument('--run-name-base', default='cnn_obs')
    parser.add_argument('--epochs', type=int, default=50, help='Number of outer epochs.')
    parser.add_argument('--inner-epochs', type=int, default=5, help='Number of inner fit epochs per outer step.')
    parser.add_argument('--batch-size', type=int, default=2000)
    parser.add_argument('--lr', type=float, default=1e-4)
    parser.add_argument('--lr-decay-factor', type=float, default=0.8)
    parser.add_argument('--log-interval', type=int, default=1)
    parser.add_argument('--end-condition', type=float, default=1e-9)
    parser.add_argument('--rounds', type=int, default=1)
    parser.add_argument('--train-ratio', type=float, default=0.7)
    parser.add_argument('--seed', type=int, default=100)
    parser.add_argument('--reg', type=float, default=0.1)
    parser.add_argument('--n-psi-train', type=int, default=100)
    parser.add_argument('--conv-filters', type=int, nargs='+', default=[32, 64, 128])
    parser.add_argument('--kernel-size', type=int, default=3)
    parser.add_argument('--pool-size', type=int, default=2)
    parser.add_argument('--dense-units', type=int, nargs='+', default=[128])
    parser.add_argument('--chunk-size', type=int, default=5000)
    parser.add_argument('--residual-form', default='projected_kv', choices=['projected_kv', 'projected_vlambda'])
    parser.add_argument('--export-checkpoint-mode', default='best', choices=['best', 'final', 'current'])
    parser.add_argument('--limit-timepoints', type=int, default=None, help='Optional debug limit for the number of time points loaded.')
    parser.add_argument('--skip-export', action='store_true')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--fresh-checkpoints', action='store_true')
    args = parser.parse_args()
    if args.output_parent is None:
        args.output_parent = resolve_local_results_root() / sanitize_tag(args.dataset_stem) / 'cnn'
    return args


ARGS = parse_args()

if ARGS.device == 'cpu':
    os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import h5py
import numpy as np
import scipy.io
import tensorflow as tf

from cnn_reskoopnet_data import (
    build_valid_transition_mask,
    load_obs_info_csv,
    make_border_aware_transition_pairs,
    reshape_obs_to_spec_and_aux,
)
from solver_resdmd_cnn import KoopmanAuxCNN, KoopmanSolver


def set_device(device_choice):
    if device_choice.lower() == 'gpu':
        gpus = tf.config.experimental.list_physical_devices('GPU')
        if gpus:
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
            print(f'Using GPU ({len(gpus)} visible devices)')
        else:
            print('GPU requested, but none are available. Falling back to CPU.')
    else:
        cpus = tf.config.experimental.list_physical_devices('CPU')
        tf.config.set_visible_devices([], 'GPU')
        tf.config.set_visible_devices(cpus, 'CPU')
        print('Using CPU')


def remove_checkpoint(checkpoint_root):
    if checkpoint_root.exists():
        shutil.rmtree(checkpoint_root)
        print(f'Removed checkpoint directory: {checkpoint_root}')
    else:
        print(f'No checkpoint directory found at: {checkpoint_root}')


def slice_data(data, start, end):
    if isinstance(data, dict):
        return {key: slice_data(value, start, end) for key, value in data.items()}
    if isinstance(data, list):
        return [slice_data(value, start, end) for value in data]
    if isinstance(data, tuple):
        return tuple(slice_data(value, start, end) for value in data)
    return data[start:end]


def index_data(data, indices):
    if isinstance(data, dict):
        return {key: index_data(value, indices) for key, value in data.items()}
    if isinstance(data, list):
        return [index_data(value, indices) for value in data]
    if isinstance(data, tuple):
        return tuple(index_data(value, indices) for value in data)
    return data[indices]


def load_obs_and_border_idx(mat_path, limit_timepoints=None):
    with h5py.File(mat_path, 'r') as f:
        obs_dataset = f['/obs']
        n_time_total = int(obs_dataset.shape[1])
        n_time = n_time_total if limit_timepoints is None else min(limit_timepoints, n_time_total)
        obs = np.array(obs_dataset[:, :n_time]).T.astype(np.float32, copy=False)
        border_idx = np.array(f['/border_idx']).reshape(-1).astype(int, copy=False)
        border_idx = border_idx[border_idx < n_time]
    return obs, border_idx, n_time_total


def restore_solver_checkpoint_for_export(solver, checkpoint_prefix, label):
    if not checkpoint_prefix:
        print(f'No {label} checkpoint path is available; keeping the current in-memory solver state.')
        return False

    checkpoint_dir = os.path.dirname(checkpoint_prefix)
    latest_checkpoint = tf.train.latest_checkpoint(checkpoint_dir)
    if latest_checkpoint is None:
        print(f'No checkpoint found under {checkpoint_dir}; keeping the current in-memory solver state.')
        return False

    K_var = tf.Variable(solver.K, trainable=False, name='K')
    reg_var = tf.Variable(solver.reg, trainable=False, name='reg')
    checkpoint = tf.train.Checkpoint(
        model=solver.model,
        K=K_var,
        eigenvectors=solver.eigenvectors_var,
        eigenvalues=solver.eigenvalues_var,
        reg=reg_var
    )
    checkpoint.restore(latest_checkpoint)

    solver.K = tf.identity(K_var.read_value())
    solver.eigenvectors = solver.eigenvectors_var.numpy().astype(np.complex128, copy=False)
    solver.eigenvalues = solver.eigenvalues_var.numpy().astype(np.complex128, copy=False)
    solver.reg = float(reg_var.numpy())
    solver.eigvec_cond = float(np.linalg.cond(solver.eigenvectors))
    solver.compute_mode()
    solver._sync_training_spectral_state()
    print(f'Restored {label} checkpoint for export: {latest_checkpoint}')
    return True


def _as_scalar_or_nan(value):
    try:
        if value is None:
            return np.nan
        arr = np.asarray(value).squeeze()
        if arr.size == 0:
            return np.nan
        if np.iscomplexobj(arr):
            arr = np.real(arr)
        return float(arr.item())
    except Exception:
        return np.nan


def _history_array(history, key):
    return np.asarray(
        [_as_scalar_or_nan(item.get(key, np.nan)) for item in history],
        dtype=np.float32
    )


def export_outputs(
        solver,
        full_inputs,
        border_idx,
        metadata,
        losses,
        val_losses,
        output_root,
        out_base,
        chunk_size,
        residual_form):
    output_root.mkdir(parents=True, exist_ok=True)
    num_chunks = int(np.ceil(full_inputs['aux'].shape[0] / float(chunk_size)))

    evalues = np.asarray(solver.eigenvalues.T, dtype=np.complex64)
    kpm_modes = np.asarray(solver.compute_mode().T, dtype=np.complex64)
    loss_array = np.asarray(losses, dtype=np.float32)
    val_loss_array = np.asarray(val_losses, dtype=np.float32)
    N_dict = int(np.shape(evalues)[0])
    outer_history = getattr(solver, 'outer_history', None) or []
    spec_group_labels = np.asarray(
        [f'{region}_{part}' for region, part in metadata['spec_group_order']],
        dtype=object,
    )

    if outer_history:
        outer_epoch_history = _history_array(outer_history, 'outer_epoch')
        outer_train_metric_history = _history_array(outer_history, 'train_metric')
        outer_val_metric_history = _history_array(outer_history, 'val_metric')
        outer_eigvec_cond_history = _history_array(outer_history, 'eigvec_cond')
        outer_lr_history = _history_array(outer_history, 'lr')
        outer_reg_history = _history_array(outer_history, 'reg')
        outer_inner_train_last_history = _history_array(outer_history, 'inner_train_last')
        outer_inner_val_last_history = _history_array(outer_history, 'inner_val_last')
    else:
        empty_history = np.array([], dtype=np.float32)
        outer_epoch_history = empty_history
        outer_train_metric_history = empty_history
        outer_val_metric_history = empty_history
        outer_eigvec_cond_history = empty_history
        outer_lr_history = empty_history
        outer_reg_history = empty_history
        outer_inner_train_last_history = empty_history
        outer_inner_val_last_history = empty_history

    shared_outputs = {
        'evalues': evalues,
        'kpm_modes': kpm_modes,
        'N_dict': N_dict,
        'loss': loss_array,
        'val_loss': val_loss_array,
        'residual_form': residual_form,
        'num_outer_epochs': int(len(outer_history)),
        'outer_epoch_history': outer_epoch_history,
        'outer_train_metric_history': outer_train_metric_history,
        'outer_val_metric_history': outer_val_metric_history,
        'outer_eigvec_cond_history': outer_eigvec_cond_history,
        'outer_lr_history': outer_lr_history,
        'outer_reg_history': outer_reg_history,
        'outer_inner_train_last_history': outer_inner_train_last_history,
        'outer_inner_val_last_history': outer_inner_val_last_history,
        'final_eigvec_cond': np.float32(_as_scalar_or_nan(getattr(solver, 'eigvec_cond', np.nan))),
        'aux_dim': np.int32(metadata['aux_dim']),
        'spec_n_freq_features': np.int32(metadata['spec_n_freq_features']),
        'spec_n_channels': np.int32(metadata['spec_n_channels']),
        'n_pairs_removed': np.int32(metadata['n_pairs_removed']),
        'n_pairs_valid': np.int32(metadata['n_pairs_valid']),
        'n_train': np.int32(metadata['n_train']),
        'n_valid': np.int32(metadata['n_valid']),
        'spec_group_labels': spec_group_labels,
        'best_checkpoint_path': str(getattr(solver, 'best_checkpoint_path', '') or ''),
        'final_checkpoint_path': str(getattr(solver, 'final_checkpoint_path', '') or ''),
    }

    summary_file = output_root / f'{out_base}_Python_resdmd_CNN_Ndict_{N_dict}_summary.mat'
    scipy.io.savemat(summary_file, {'EDMD_outputs': shared_outputs})
    print(f'saved summary: {summary_file}')

    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, full_inputs['aux'].shape[0])
        chunk_inputs = slice_data(full_inputs, start_idx, end_idx)
        efuns = np.asarray(solver.eigenfunctions(chunk_inputs), dtype=np.complex64)

        out_file = output_root / f'{out_base}_Python_resdmd_CNN_Ndict_{N_dict}_outputs_{i + 1}.mat'
        scipy.io.savemat(out_file, {'EDMD_outputs': {'efuns': efuns, **shared_outputs}})
        print(f'saved: {out_file}')

        del chunk_inputs
        del efuns
        gc.collect()

    all_valid_pair_idx = np.flatnonzero(
        build_valid_transition_mask(full_inputs['aux'].shape[0], border_idx=border_idx)
    )
    num_pair_chunks = int(np.ceil(all_valid_pair_idx.size / float(chunk_size)))

    for i in range(num_pair_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, all_valid_pair_idx.size)
        idx_chunk = all_valid_pair_idx[start_idx:end_idx]
        x_chunk = index_data(full_inputs, idx_chunk)
        y_chunk = index_data(full_inputs, idx_chunk + 1)

        psi_x = np.asarray(solver.dic.call(x_chunk).numpy().T, dtype=np.float32)
        psi_y = np.asarray(solver.dic.call(y_chunk).numpy().T, dtype=np.float32)

        out_file = output_root / f'{out_base}_Python_resdmd_CNN_Ndict_{N_dict}_outputs_Psi_{i + 1}.mat'
        scipy.io.savemat(out_file, {'EDMD_outputs': {'Psi_X': psi_x, 'Psi_Y': psi_y}})
        print(f'saved: {out_file}')

        del idx_chunk
        del x_chunk
        del y_chunk
        del psi_x
        del psi_y
        gc.collect()


def main():
    set_device(ARGS.device)

    mat_filename = ARGS.data_filename or f'{ARGS.dataset_stem}_low50_high250_g2_{ARGS.observable_mode}_single.mat'
    info_csv_filename = ARGS.info_csv_filename or f'{ARGS.dataset_stem}_low50_high250_g2_{ARGS.observable_mode}_single_obs_info.csv'
    mat_path = ARGS.data_root / mat_filename
    info_csv_path = ARGS.data_root / info_csv_filename

    run_label = f'{ARGS.run_name_base}_{ARGS.residual_form}_{ARGS.observable_mode}_{ARGS.device}'
    run_root = ARGS.output_parent / run_label
    checkpoint_root = run_root / 'checkpoints'
    output_root = run_root / 'outputs'
    log_root = run_root / 'logs'
    checkpoint_path = checkpoint_root / 'final' / 'model.ckpt'
    best_checkpoint_path = checkpoint_root / 'best' / 'model.ckpt'

    print(f'mat_path: {mat_path}')
    print(f'info_csv_path: {info_csv_path}')
    print(f'run_root: {run_root}')

    if ARGS.fresh_checkpoints:
        remove_checkpoint(checkpoint_root)

    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    best_checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    output_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)

    obs, border_idx, n_time_total = load_obs_and_border_idx(mat_path, ARGS.limit_timepoints)
    obs_info_rows = load_obs_info_csv(info_csv_path)

    print(f'Loaded obs shape: {obs.shape} (source total time points: {n_time_total})')
    print(f'Loaded border_idx count: {border_idx.size}')

    full_inputs, reshape_meta = reshape_obs_to_spec_and_aux(
        obs,
        obs_info_rows,
        output_structure='dict'
    )
    data_train, data_valid, pair_meta = make_border_aware_transition_pairs(
        full_inputs,
        border_idx=border_idx,
        train_ratio=ARGS.train_ratio,
        seed=ARGS.seed
    )

    metadata = {**reshape_meta, **pair_meta}
    print(f'Dataset metadata: {metadata}')
    print(f"train spec shape: {data_train[0]['spec'].shape}, aux shape: {data_train[0]['aux'].shape}")
    print(f"valid spec shape: {data_valid[0]['spec'].shape}, aux shape: {data_valid[0]['aux'].shape}")

    basis_function = KoopmanAuxCNN(
        spec_input_shape=(metadata['spec_n_freq_features'], metadata['spec_n_channels']),
        aux_dim=metadata['aux_dim'],
        conv_filters=tuple(ARGS.conv_filters),
        kernel_size=ARGS.kernel_size,
        pool_size=ARGS.pool_size,
        dense_units=tuple(ARGS.dense_units),
        n_psi_train=ARGS.n_psi_train,
    )

    solver = KoopmanSolver(
        dic=basis_function,
        target_dim=metadata['aux_dim'],
        reg=ARGS.reg,
        residual_form=ARGS.residual_form,
        training_policy='float32' if ARGS.device == 'cpu' else None,
    )

    loss_all = []
    val_loss_all = []
    completed_rounds = 0
    stop_flag = False

    for n_round in range(ARGS.rounds):
        print(f'Round number: {n_round}')
        temp_loss, temp_val_loss, stop_flag, _ = solver.build(
            data_train=data_train,
            data_valid=data_valid,
            epochs=ARGS.epochs,
            batch_size=ARGS.batch_size,
            lr=ARGS.lr,
            log_interval=ARGS.log_interval,
            lr_decay_factor=ARGS.lr_decay_factor,
            Nepoch=ARGS.inner_epochs,
            end_condition=ARGS.end_condition,
            checkpoint_path=str(checkpoint_path),
            best_checkpoint_path=str(best_checkpoint_path),
            resume=ARGS.resume,
        )
        loss_all.extend(temp_loss)
        val_loss_all.extend(temp_val_loss)
        completed_rounds = n_round + 1
        if stop_flag:
            break

    print(f'Completed rounds: {completed_rounds}, stop_flag={stop_flag}')
    print(f'Collected losses: train={len(loss_all)}, valid={len(val_loss_all)}')

    if ARGS.skip_export:
        return

    if ARGS.export_checkpoint_mode == 'best':
        restore_solver_checkpoint_for_export(solver, str(best_checkpoint_path), 'best')
    elif ARGS.export_checkpoint_mode == 'final':
        restore_solver_checkpoint_for_export(solver, str(checkpoint_path), 'final')
    else:
        print('Keeping current in-memory solver state for export.')

    export_outputs(
        solver=solver,
        full_inputs=full_inputs,
        border_idx=border_idx,
        metadata=metadata,
        losses=loss_all,
        val_losses=val_loss_all,
        output_root=output_root,
        out_base=mat_path.stem,
        chunk_size=ARGS.chunk_size,
        residual_form=ARGS.residual_form,
    )


if __name__ == '__main__':
    main()
