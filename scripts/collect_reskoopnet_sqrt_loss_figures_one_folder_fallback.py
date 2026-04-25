from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
import numpy as np
from scipy.io import loadmat

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SOURCE_ROOT = Path("/mnt/e/autodl_results")
REPO_ROOT = Path("/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished")
OUTPUT_ROOT = REPO_ROOT / "analysis_outputs" / "reskoopnet_training_sqrt_loss_figures_one_folder_20260424"


def to_text(value, default=""):
    if value is None:
        return str(default)
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore")
    if isinstance(value, str):
        return value
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return str(default)
        return to_text(value.reshape(-1)[0], default)
    return str(value)


def get_field(obj, name, default=None):
    return getattr(obj, name, default) if hasattr(obj, name) else default


def num_array(obj, name):
    value = get_field(obj, name, None)
    if value is None:
        return np.array([], dtype=float)
    arr = np.asarray(value)
    if arr.size == 0:
        return np.array([], dtype=float)
    try:
        return arr.astype(float).ravel()
    except Exception:
        return np.array([], dtype=float)


def first_num_array(obj, names):
    for name in names:
        arr = num_array(obj, name)
        if arr.size:
            return arr
    return np.array([], dtype=float)


def sqrt_nonneg(arr):
    arr = np.asarray(arr, dtype=float).copy().ravel()
    if arr.size == 0:
        return arr
    finite = np.isfinite(arr)
    arr[finite] = np.sqrt(np.maximum(arr[finite], 0.0))
    return arr


def family_for(path_text: str) -> str:
    return "BOLD" if "/bold/" in path_text.lower() else "BLP"


def dataset_for(path_text: str) -> str:
    lower = path_text.lower()
    for token in ("e10gb1", "e10fv1", "e10gh1", "f12m01"):
        if f"/{token}/" in lower:
            return token
    return "unknown"


def observable_for(path_text: str) -> str:
    lower = path_text.lower()
    if "complex_split" in lower:
        return "complex_split"
    if "/bold/" in lower:
        return "bold"
    if "abs" in lower:
        return "abs"
    return "unknown"


def residual_for(path_text: str) -> str:
    lower = path_text.lower()
    if "projected_vlambda" in lower:
        return "projected_vlambda"
    if "projected_kv" in lower:
        return "projected_kv"
    return "unknown"


def sanitize(text: str) -> str:
    out = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text)
    while "__" in out:
        out = out.replace("__", "_")
    out = out.strip("_")
    return out or "untitled"


def plot_series(ax, x, y, label):
    if len(x) == 0 or len(y) == 0:
        return
    n = min(len(x), len(y))
    ax.plot(x[:n], y[:n], linewidth=1.2, marker=None if n > 20 else "o", label=label)


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    summary_files = sorted(SOURCE_ROOT.rglob("*_summary.mat"))
    records = []
    export_count = 0
    skip_count = 0

    for summary_path in summary_files:
        path_text = summary_path.as_posix()
        family = family_for(path_text)
        dataset = dataset_for(path_text)
        observable = observable_for(path_text)
        residual = residual_for(path_text)
        run_label = summary_path.parent.name
        loss_mode = ""
        requested_loss_mode = ""
        diag_png = ""
        metrics_png = ""
        status = "skipped"
        reason = ""

        try:
            loaded = loadmat(str(summary_path), squeeze_me=True, struct_as_record=False)
            if "EDMD_outputs" not in loaded:
                raise RuntimeError("EDMD_outputs missing")
            edmd = loaded["EDMD_outputs"]
            loss_mode = to_text(get_field(edmd, "loss_mode", "squared"), "squared")
            requested_loss_mode = to_text(get_field(edmd, "requested_loss_mode", loss_mode), loss_mode)
            dataset = to_text(get_field(edmd, "dataset_stem", dataset), dataset)
            observable = to_text(get_field(edmd, "observable_mode", observable), observable)
            residual = to_text(get_field(edmd, "residual_form", residual), residual)
            run_label = Path(to_text(get_field(edmd, "run_label", run_label), run_label)).name

            if loss_mode.lower() != "squared":
                reason = "non_squared_loss_mode"
            else:
                inner_loss = num_array(edmd, "loss")
                inner_val_loss = num_array(edmd, "val_loss")
                if inner_loss.size == 0 and inner_val_loss.size == 0:
                    reason = "missing_inner_loss"
                else:
                    outer_epoch = num_array(edmd, "outer_epoch_history")
                    outer_train_sq = first_num_array(edmd, ["outer_train_metric_squared_history", "outer_train_metric_history"])
                    outer_val_sq = first_num_array(edmd, ["outer_val_metric_squared_history", "outer_val_metric_history"])
                    outer_cond = num_array(edmd, "outer_eigvec_cond_history")
                    inner_train_sq = first_num_array(edmd, ["inner_train_metric_squared_history", "loss"])
                    inner_val_sq = first_num_array(edmd, ["inner_val_metric_squared_history", "val_loss"])
                    inner_train_per_dim = num_array(edmd, "inner_train_metric_per_dim_history")
                    inner_val_per_dim = num_array(edmd, "inner_val_metric_per_dim_history")

                    if outer_epoch.size == 0:
                        outer_epoch = np.arange(1, max(outer_train_sq.size, outer_val_sq.size) + 1, dtype=float)

                    sqrt_inner_loss = sqrt_nonneg(inner_loss)
                    sqrt_inner_val_loss = sqrt_nonneg(inner_val_loss)
                    sqrt_outer_train = sqrt_nonneg(outer_train_sq)
                    sqrt_outer_val = sqrt_nonneg(outer_val_sq)
                    sqrt_inner_sq = sqrt_nonneg(inner_train_sq)
                    sqrt_inner_val_sq = sqrt_nonneg(inner_val_sq)
                    sqrt_inner_per_dim = sqrt_nonneg(inner_train_per_dim)
                    sqrt_inner_val_per_dim = sqrt_nonneg(inner_val_per_dim)

                    base_tag = sanitize("__".join([family, dataset, residual, observable, run_label]))

                    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
                    fig.suptitle(f"{dataset} {observable} {residual} sqrt squared residual diagnostics")

                    ax = axes[0, 0]
                    plot_series(ax, np.arange(1, sqrt_inner_loss.size + 1), sqrt_inner_loss, "train")
                    plot_series(ax, np.arange(1, sqrt_inner_val_loss.size + 1), sqrt_inner_val_loss, "val")
                    ax.set_xlabel("inner fit epoch")
                    ax.set_ylabel("sqrt(loss)")
                    ax.set_title("Inner RMS residual")
                    ax.grid(True, alpha=0.3)
                    ax.legend()

                    ax = axes[0, 1]
                    plot_series(ax, outer_epoch, sqrt_outer_train, "train")
                    plot_series(ax, outer_epoch, sqrt_outer_val, "val")
                    ax.set_xlabel("outer epoch")
                    ax.set_ylabel("sqrt(outer squared metric)")
                    ax.set_title("Outer RMS residual")
                    ax.grid(True, alpha=0.3)
                    ax.legend()

                    ax = axes[1, 0]
                    if outer_cond.size:
                        n = min(outer_epoch.size, outer_cond.size)
                        ax.plot(outer_epoch[:n], np.log10(np.maximum(outer_cond[:n], np.finfo(float).tiny)), marker="o", linewidth=1.2)
                    ax.set_xlabel("outer epoch")
                    ax.set_ylabel("log10(eigvec cond)")
                    ax.set_title("Spectral conditioning")
                    ax.grid(True, alpha=0.3)

                    ax = axes[1, 1]
                    if sqrt_outer_val.size:
                        n = min(outer_epoch.size, sqrt_outer_val.size)
                        epoch_scatter = outer_epoch[:n]
                        val_scatter = sqrt_outer_val[:n]
                        finite = np.isfinite(val_scatter)
                        if finite.any():
                            ax.scatter(epoch_scatter[finite], val_scatter[finite], s=24, label="outer val")
                            finite_idx = np.where(finite)[0]
                            best_local = finite_idx[np.argmin(val_scatter[finite])]
                            ax.scatter([epoch_scatter[best_local]], [val_scatter[best_local]], s=80, label=f"best epoch {int(round(epoch_scatter[best_local]))}")
                            ax.legend()
                    ax.set_xlabel("outer epoch")
                    ax.set_ylabel("sqrt(val outer squared metric)")
                    ax.set_title("Best validation epoch")
                    ax.grid(True, alpha=0.3)

                    diag_png_path = OUTPUT_ROOT / f"{base_tag}__sqrt_loss_diagnostics.png"
                    diag_pdf_path = OUTPUT_ROOT / f"{base_tag}__sqrt_loss_diagnostics.pdf"
                    fig.savefig(diag_png_path, dpi=180)
                    fig.savefig(diag_pdf_path)
                    plt.close(fig)
                    diag_png = str(diag_png_path)

                    fig, axes = plt.subplots(2, 1, figsize=(12, 7.6), constrained_layout=True)
                    fig.suptitle(f"{dataset} {observable} {residual} sqrt inner metrics")

                    ax = axes[0]
                    plot_series(ax, np.arange(1, sqrt_inner_sq.size + 1), sqrt_inner_sq, "train")
                    plot_series(ax, np.arange(1, sqrt_inner_val_sq.size + 1), sqrt_inner_val_sq, "val")
                    ax.set_xlabel("inner fit epoch")
                    ax.set_ylabel("sqrt(squared)")
                    ax.set_title("Inner RMS residual from squared metric")
                    ax.grid(True, alpha=0.3)
                    ax.legend()

                    ax = axes[1]
                    plot_series(ax, np.arange(1, sqrt_inner_per_dim.size + 1), sqrt_inner_per_dim, "train")
                    plot_series(ax, np.arange(1, sqrt_inner_val_per_dim.size + 1), sqrt_inner_val_per_dim, "val")
                    ax.set_xlabel("inner fit epoch")
                    ax.set_ylabel("sqrt(per_dim)")
                    ax.set_title("Inner per-dimension RMS residual")
                    ax.grid(True, alpha=0.3)
                    ax.legend()

                    metrics_png_path = OUTPUT_ROOT / f"{base_tag}__sqrt_inner_metrics.png"
                    metrics_pdf_path = OUTPUT_ROOT / f"{base_tag}__sqrt_inner_metrics.pdf"
                    fig.savefig(metrics_png_path, dpi=180)
                    fig.savefig(metrics_pdf_path)
                    plt.close(fig)
                    metrics_png = str(metrics_png_path)

                    status = "exported"
                    export_count += 1
        except Exception as exc:
            reason = str(exc)

        if status != "exported":
            skip_count += 1

        records.append(
            {
                "status": status,
                "reason": reason,
                "summary_path": str(summary_path),
                "family": family,
                "dataset_stem": dataset,
                "observable_mode": observable,
                "residual_form": residual,
                "run_label": run_label,
                "loss_mode": loss_mode,
                "requested_loss_mode": requested_loss_mode,
                "sqrt_diagnostics_png": diag_png,
                "sqrt_inner_metrics_png": metrics_png,
            }
        )

    fieldnames = [
        "status",
        "reason",
        "summary_path",
        "family",
        "dataset_stem",
        "observable_mode",
        "residual_form",
        "run_label",
        "loss_mode",
        "requested_loss_mode",
        "sqrt_diagnostics_png",
        "sqrt_inner_metrics_png",
    ]
    with (OUTPUT_ROOT / "sqrt_loss_figures_index.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    (OUTPUT_ROOT / "README.txt").write_text(
        "ResKoopNet sqrt loss figures (one folder)\n"
        f"Source root: {SOURCE_ROOT}\n"
        f"Exported runs: {export_count}\n"
        f"Skipped runs: {skip_count}\n"
        "Included loss mode: squared only\n",
        encoding="utf-8",
    )

    print(f"Saved sqrt loss figures to: {OUTPUT_ROOT}")
    print(f"Exported: {export_count} | Skipped: {skip_count}")


if __name__ == "__main__":
    main()
