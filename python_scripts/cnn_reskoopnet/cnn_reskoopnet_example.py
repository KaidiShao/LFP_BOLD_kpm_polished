import numpy as np
from solver_resdmd_cnn import KoopmanCNN, KoopmanSolver

# Example 1: spectrogram snapshots with shape [num_snapshots, freq_bins, time_bins, channels].
# Add a singleton channel axis if your data is [N, H, W].
spectrogram_data = np.random.randn(1000, 64, 64, 1).astype("float32")

data_train = [spectrogram_data[:-1], spectrogram_data[1:]]
data_valid = [spectrogram_data[:-1], spectrogram_data[1:]]

dic = KoopmanCNN(
    input_shape=(64, 64, 1),
    conv_filters=(16, 32, 64),
    kernel_size=3,
    pool_size=2,
    dense_units=(128,),
    n_psi_train=32,
    include_raw_state=True,
)

solver = KoopmanSolver(
    dic=dic,
    target_dim=(64, 64, 1),
    reg=1e-3,
)

losses, val_losses, stop_flag, learning_rate_changes = solver.build(
    data_train=data_train,
    data_valid=data_valid,
    epochs=10,
    batch_size=32,
    lr=1e-3,
    log_interval=1,
    lr_decay_factor=0.5,
    Nepoch=5,
)

# Example 2: if your BOLD data is a 3D volume per snapshot, use input_shape=(D, H, W, C).
# The same KoopmanCNN class will switch to Conv3D automatically as long as the channel axis is present.
