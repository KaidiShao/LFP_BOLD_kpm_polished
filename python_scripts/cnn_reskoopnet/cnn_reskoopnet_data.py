import csv
from collections import OrderedDict

import numpy as np


def load_obs_info_csv(csv_path):
    with open(csv_path, newline='') as f:
        rows = list(csv.DictReader(f))

    rows.sort(key=lambda row: int(row['dim_id']))
    return rows


def _leading_dim(data):
    if isinstance(data, dict):
        first_leaf = next(iter(data.values()))
        return _leading_dim(first_leaf)
    if isinstance(data, (list, tuple)):
        return _leading_dim(data[0])
    return int(np.shape(data)[0])


def _index_data(data, indices):
    if isinstance(data, dict):
        return {key: _index_data(value, indices) for key, value in data.items()}
    if isinstance(data, list):
        return [_index_data(value, indices) for value in data]
    if isinstance(data, tuple):
        return tuple(_index_data(value, indices) for value in data)
    return data[indices]


def build_valid_transition_mask(length, border_idx=None, one_based=True):
    if length < 2:
        return np.zeros((0,), dtype=bool)

    mask = np.ones((length - 1,), dtype=bool)
    if border_idx is None:
        return mask

    border_idx = np.asarray(border_idx).reshape(-1)
    if border_idx.size == 0:
        return mask

    border_idx = border_idx.astype(int, copy=False)
    invalid_pair_idx = border_idx - 1 if one_based else border_idx
    invalid_pair_idx = invalid_pair_idx[
        (invalid_pair_idx >= 0) & (invalid_pair_idx < (length - 1))
    ]
    mask[invalid_pair_idx] = False
    return mask


def reshape_obs_to_spec_and_aux(
        obs,
        obs_info_rows,
        output_structure='dict',
        extra_aux=None):
    obs = np.asarray(obs)
    if obs.ndim != 2:
        raise ValueError("obs must be a 2D array shaped [time, features].")

    if obs.shape[1] != len(obs_info_rows):
        raise ValueError(
            f"obs has {obs.shape[1]} columns, but obs_info has {len(obs_info_rows)} rows."
        )

    raw_indices = []
    spec_groups = OrderedDict()

    for col_idx, row in enumerate(obs_info_rows):
        if row['source'] == 'blp':
            raw_indices.append(col_idx)
            continue

        if row['source'] != 'spectrogram':
            continue

        group_key = (row['region_label'], row['part'])
        spec_groups.setdefault(group_key, []).append(col_idx)

    if not raw_indices:
        raise ValueError("No raw BLP observables were found in obs_info.")
    if not spec_groups:
        raise ValueError("No spectrogram observables were found in obs_info.")

    aux = obs[:, raw_indices].astype(np.float32, copy=False)

    spec_blocks = []
    for _, cols in spec_groups.items():
        spec_blocks.append(obs[:, cols].astype(np.float32, copy=False))

    spec_block_lengths = [block.shape[1] for block in spec_blocks]
    if len(set(spec_block_lengths)) != 1:
        raise ValueError(
            "Spectrogram groups do not share the same frequency dimension. "
            f"Found lengths: {spec_block_lengths}"
        )

    spec = np.stack(spec_blocks, axis=-1)  # [time, freq_features, spec_channels]

    if extra_aux is not None:
        extra_aux = np.asarray(extra_aux, dtype=np.float32)
        if extra_aux.ndim == 1:
            extra_aux = extra_aux[:, None]
        if extra_aux.shape[0] != aux.shape[0]:
            raise ValueError(
                f"extra_aux has {extra_aux.shape[0]} rows, expected {aux.shape[0]}."
            )
        aux = np.concatenate([aux, extra_aux], axis=-1)

    metadata = {
        'raw_indices': raw_indices,
        'spec_group_order': list(spec_groups.keys()),
        'spec_n_freq_features': int(spec.shape[1]),
        'spec_n_channels': int(spec.shape[2]),
        'aux_dim': int(aux.shape[1])
    }

    if output_structure == 'dict':
        inputs = {'spec': spec, 'aux': aux}
    elif output_structure == 'tuple':
        inputs = (spec, aux)
    else:
        raise ValueError("output_structure must be either 'dict' or 'tuple'.")

    return inputs, metadata


def make_border_aware_transition_pairs(
        inputs,
        border_idx=None,
        train_ratio=0.7,
        seed=100,
        shuffle=True,
        one_based_border_idx=True):
    n_time = _leading_dim(inputs)
    valid_mask = build_valid_transition_mask(
        n_time,
        border_idx=border_idx,
        one_based=one_based_border_idx
    )
    pair_indices = np.flatnonzero(valid_mask)

    rng = np.random.RandomState(seed)
    if shuffle:
        pair_indices = rng.permutation(pair_indices)

    split_idx = int(train_ratio * pair_indices.size)
    train_indices = pair_indices[:split_idx]
    valid_indices = pair_indices[split_idx:]

    data_train = [
        _index_data(inputs, train_indices),
        _index_data(inputs, train_indices + 1)
    ]
    data_valid = [
        _index_data(inputs, valid_indices),
        _index_data(inputs, valid_indices + 1)
    ]

    metadata = {
        'n_time': int(n_time),
        'n_pairs_total': int(max(n_time - 1, 0)),
        'n_pairs_valid': int(pair_indices.size),
        'n_pairs_removed': int((n_time - 1) - pair_indices.size),
        'n_train': int(train_indices.size),
        'n_valid': int(valid_indices.size),
        'train_indices': train_indices,
        'valid_indices': valid_indices
    }
    return data_train, data_valid, metadata


def build_obs_1dcnn_dataset(
        obs,
        obs_info_rows,
        border_idx=None,
        train_ratio=0.7,
        seed=100,
        output_structure='dict',
        extra_aux=None):
    inputs, reshape_meta = reshape_obs_to_spec_and_aux(
        obs,
        obs_info_rows,
        output_structure=output_structure,
        extra_aux=extra_aux
    )
    data_train, data_valid, pair_meta = make_border_aware_transition_pairs(
        inputs,
        border_idx=border_idx,
        train_ratio=train_ratio,
        seed=seed
    )

    metadata = {}
    metadata.update(reshape_meta)
    metadata.update(pair_meta)
    return data_train, data_valid, metadata
