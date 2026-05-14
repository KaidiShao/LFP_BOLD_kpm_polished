from __future__ import annotations

import json
import math
import sys
import time
from datetime import datetime
from pathlib import Path
from xml.sax.saxutils import escape


RESULT_BASE = Path("E:/autodl_results_local/bold_wsl")
REPORT_DIR = Path("tmp/monitor_reports")
STALE_THRESHOLD = 2000

TASKS = [
    ("E10.gH1", "e10gh1", "gsvd100_ds", "bold_batch4_gsvd100_ds_inner1_lr1e5_test1000_20260513_e10gh1"),
    ("E10.gb1", "e10gb1", "global_svd100", "bold_batch4_global_svd100_inner1_lr1e5_reg0p1_untilstale_20260513_e10gb1"),
    ("E10.gb1", "e10gb1", "gsvd100_ds", "bold_batch4_gsvd100_ds_inner1_lr1e5_reg0p1_untilstale_20260513_e10gb1"),
    ("E10.fV1", "e10fV1", "global_svd100", "bold_batch4_global_svd100_inner1_lr1e5_reg0p1_untilstale_20260513_e10fV1"),
    ("E10.fV1", "e10fV1", "gsvd100_ds", "bold_batch4_gsvd100_ds_inner1_lr1e5_reg0p1_untilstale_20260513_e10fV1"),
    ("F12.m01", "f12m01", "global_svd100", "bold_batch4_global_svd100_inner1_lr1e5_reg0p1_untilstale_20260513_f12m01"),
    ("F12.m01", "f12m01", "gsvd100_ds", "bold_batch4_gsvd100_ds_inner1_lr1e5_reg0p1_untilstale_20260513_f12m01"),
    ("E10.gW1", "e10gw1", "global_svd100", "bold_batch4_global_svd100_inner1_lr1e5_reg0p1_untilstale_20260513_e10gw1"),
    ("E10.gW1", "e10gw1", "gsvd100_ds", "bold_batch4_gsvd100_ds_inner1_lr1e5_reg0p1_untilstale_20260513_e10gw1"),
]


def run_label(experiment: str, observable: str) -> str:
    return f"mlp_obs_{experiment}_projected_vlambda_{observable}"


def state_path(stem: str, experiment: str, observable: str) -> Path:
    return (
        RESULT_BASE
        / stem
        / "mlp"
        / "checkpoints"
        / run_label(experiment, observable)
        / "final"
        / "training_state.json"
    )


def read_json_retry(path: Path) -> dict:
    last_error: Exception | None = None
    for _ in range(8):
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as exc:
            last_error = exc
            time.sleep(0.25)
    if last_error:
        raise last_error
    raise RuntimeError(f"failed to read {path}")


def summarize_task(dataset: str, stem: str, observable: str, experiment: str) -> tuple[dict, dict | None]:
    path = state_path(stem, experiment, observable)
    base = {
        "dataset": dataset,
        "stem": stem,
        "observable": observable,
        "experiment": experiment,
        "label": run_label(experiment, observable),
        "path": str(path),
    }
    if not path.exists():
        return {**base, "state": "not_started"}, None

    try:
        state = read_json_retry(path)
    except Exception as exc:
        return {**base, "state": "read_error", "error": str(exc)}, None

    history = state.get("outer_history") or []
    if not history:
        return {**base, "state": "empty"}, None

    epochs = [int(row["outer_epoch"]) for row in history]
    train = [float(row["train_metric"]) for row in history]
    val = [float(row["val_metric"]) for row in history]
    best_idx = min(range(len(val)), key=val.__getitem__)
    best_epoch = epochs[best_idx]
    final_epoch = epochs[-1]
    stale_epochs = final_epoch - best_epoch
    status = "done" if stale_epochs >= STALE_THRESHOLD else "running_or_recent"

    best_so_far = []
    current = math.inf
    for value in val:
        current = min(current, value)
        best_so_far.append(current)

    summary = {
        **base,
        "state": status,
        "epochs": len(history),
        "best_epoch": best_epoch,
        "best_val": val[best_idx],
        "best_train": train[best_idx],
        "final_epoch": final_epoch,
        "final_val": val[-1],
        "final_train": train[-1],
        "stale_epochs": stale_epochs,
    }
    series = {"epoch": epochs, "val": val, "best": best_so_far, "best_idx": best_idx}
    return summary, series


def format_val(value: object) -> str:
    if not isinstance(value, (float, int)):
        return "-"
    if abs(float(value)) >= 1000:
        return f"{float(value):.3g}"
    return f"{float(value):.4g}"


def compare_changes(previous: list[dict], current: list[dict]) -> list[str]:
    previous_by_key = {(r.get("dataset"), r.get("observable")): r for r in previous}
    changes = []
    for row in current:
        key = (row.get("dataset"), row.get("observable"))
        old = previous_by_key.get(key)
        name = f"{row['dataset']} {row['observable']}"
        if old is None:
            if row["state"] != "not_started":
                changes.append(f"- {name}: started at epoch {row.get('final_epoch', '-')}.")
            continue
        if old.get("state") != row.get("state"):
            changes.append(f"- {name}: state {old.get('state')} -> {row.get('state')}.")
        if row.get("best_epoch") and old.get("best_epoch") != row.get("best_epoch"):
            changes.append(
                f"- {name}: new best {format_val(row.get('best_val'))} @ {row.get('best_epoch')}."
            )
        elif row.get("final_epoch") and old.get("final_epoch") != row.get("final_epoch"):
            changes.append(
                f"- {name}: advanced to epoch {row.get('final_epoch')} "
                f"(current {format_val(row.get('final_val'))}, stale {row.get('stale_epochs')})."
            )
    return changes


def polyline(points: list[tuple[float, float]]) -> str:
    return " ".join(f"{x:.1f},{y:.1f}" for x, y in points)


def downsample(values: list[tuple[float, float]], max_points: int = 700) -> list[tuple[float, float]]:
    if len(values) <= max_points:
        return values
    step = max(1, math.ceil(len(values) / max_points))
    sampled = values[::step]
    if sampled[-1] != values[-1]:
        sampled.append(values[-1])
    return sampled


def render_svg(series_items: list[tuple[dict, dict]], path: Path) -> None:
    if not series_items:
        path.write_text("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"120\" />", encoding="utf-8")
        return

    cols = 2
    panel_w = 560
    panel_h = 245
    margin = 52
    title_h = 42
    rows = math.ceil(len(series_items) / cols)
    width = cols * panel_w
    height = rows * panel_h + 56
    parts = [
        f"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\">",
        "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>",
        "<style>text{font-family:Arial,sans-serif;fill:#222}.small{font-size:12px}.title{font-size:13px;font-weight:700}.axis{stroke:#999;stroke-width:1}.grid{stroke:#ddd;stroke-width:1}.val{fill:none;stroke:#4C78A8;stroke-width:1.2;opacity:.85}.best{fill:none;stroke:#1B7F3A;stroke-width:2.0}.dot{fill:#D62728}</style>",
        "<text x=\"16\" y=\"24\" class=\"title\">BOLD global_svd100 / gsvd100_ds monitor</text>",
    ]

    for idx, (summary, series) in enumerate(series_items):
        col = idx % cols
        row = idx // cols
        x0 = col * panel_w
        y0 = 42 + row * panel_h
        plot_x = x0 + margin
        plot_y = y0 + title_h + 8
        plot_w = panel_w - margin - 24
        plot_h = panel_h - title_h - 44

        epochs = series["epoch"]
        vals = [max(float(v), 1e-12) for v in series["val"]]
        bests = [max(float(v), 1e-12) for v in series["best"]]
        xs = epochs
        y_logs = [math.log10(v) for v in vals + bests if math.isfinite(v) and v > 0]
        if not y_logs:
            continue
        xmin, xmax = min(xs), max(xs)
        if xmin == xmax:
            xmax = xmin + 1
        ymin, ymax = min(y_logs), max(y_logs)
        if abs(ymax - ymin) < 1e-9:
            ymax = ymin + 1.0

        def map_x(epoch: float) -> float:
            return plot_x + (epoch - xmin) / (xmax - xmin) * plot_w

        def map_y(value: float) -> float:
            logv = math.log10(max(value, 1e-12))
            return plot_y + (ymax - logv) / (ymax - ymin) * plot_h

        val_points = downsample([(map_x(e), map_y(v)) for e, v in zip(epochs, vals)])
        best_points = downsample([(map_x(e), map_y(v)) for e, v in zip(epochs, bests)])
        best_idx = int(series["best_idx"])
        bx, by = map_x(epochs[best_idx]), map_y(vals[best_idx])

        title = (
            f"{summary['dataset']} {summary['observable']} | "
            f"best {format_val(summary['best_val'])} @ {summary['best_epoch']} | "
            f"final {format_val(summary['final_val'])} @ {summary['final_epoch']} | "
            f"stale {summary['stale_epochs']}"
        )
        parts.extend(
            [
                f"<text x=\"{x0 + 12}\" y=\"{y0 + 18}\" class=\"title\">{escape(title)}</text>",
                f"<rect x=\"{plot_x}\" y=\"{plot_y}\" width=\"{plot_w}\" height=\"{plot_h}\" fill=\"#fafafa\" stroke=\"#ccc\"/>",
                f"<line x1=\"{plot_x}\" y1=\"{plot_y + plot_h / 2}\" x2=\"{plot_x + plot_w}\" y2=\"{plot_y + plot_h / 2}\" class=\"grid\"/>",
                f"<polyline class=\"val\" points=\"{polyline(val_points)}\"/>",
                f"<polyline class=\"best\" points=\"{polyline(best_points)}\"/>",
                f"<circle class=\"dot\" cx=\"{bx:.1f}\" cy=\"{by:.1f}\" r=\"3.2\"/>",
                f"<text x=\"{plot_x}\" y=\"{plot_y + plot_h + 18}\" class=\"small\">epoch {xmin} - {xmax}</text>",
                f"<text x=\"{plot_x + plot_w - 110}\" y=\"{plot_y + plot_h + 18}\" class=\"small\">log val loss</text>",
            ]
        )

    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def make_markdown(rows: list[dict], changes: list[str], svg_path: Path, generated_at: str) -> str:
    active = [r for r in rows if r["state"] == "running_or_recent"]
    done = [r for r in rows if r["state"] == "done"]
    pending = [r for r in rows if r["state"] == "not_started"]

    lines = [
        f"**BOLD Monitor {generated_at}**",
        "",
    ]
    if changes:
        lines.append("Changes since last check:")
        lines.extend(changes[:8])
        lines.append("")
    else:
        lines.append("No major state change since last check.")
        lines.append("")

    if active:
        lines.append("Active / not yet stale:")
        for row in active:
            lines.append(
                f"- {row['dataset']} {row['observable']}: epoch {row['final_epoch']}, "
                f"best {format_val(row['best_val'])} @ {row['best_epoch']}, "
                f"current {format_val(row['final_val'])}, stale {row['stale_epochs']}"
            )
        lines.append("")

    if done:
        lines.append("Done by 2000-stale rule:")
        for row in done:
            lines.append(
                f"- {row['dataset']} {row['observable']}: best {format_val(row['best_val'])} "
                f"@ {row['best_epoch']}, final {format_val(row['final_val'])} @ {row['final_epoch']}"
            )
        lines.append("")

    if pending:
        lines.append("Not started:")
        lines.append("- " + ", ".join(f"{r['dataset']} {r['observable']}" for r in pending))
        lines.append("")

    abs_svg = svg_path.resolve().as_posix()
    lines.append(f"![latest loss curves]({abs_svg})")
    return "\n".join(lines)


def main() -> int:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    latest_json = REPORT_DIR / "bold_global_svd_monitor_latest.json"
    previous = []
    if latest_json.exists():
        try:
            previous = json.loads(latest_json.read_text(encoding="utf-8"))
        except Exception:
            previous = []

    rows = []
    series_items = []
    for task in TASKS:
        summary, series = summarize_task(*task)
        rows.append(summary)
        if series is not None:
            series_items.append((summary, series))

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    svg_path = REPORT_DIR / f"bold_global_svd_loss_curves_{stamp}.svg"
    latest_svg = REPORT_DIR / "bold_global_svd_loss_curves_latest.svg"
    render_svg(series_items, svg_path)
    latest_svg.write_text(svg_path.read_text(encoding="utf-8"), encoding="utf-8")

    changes = compare_changes(previous, rows)
    markdown = make_markdown(rows, changes, svg_path, generated_at)
    md_path = REPORT_DIR / f"bold_global_svd_monitor_{stamp}.md"
    latest_md = REPORT_DIR / "bold_global_svd_monitor_latest.md"
    md_path.write_text(markdown, encoding="utf-8")
    latest_md.write_text(markdown, encoding="utf-8")
    latest_json.write_text(json.dumps(rows, indent=2, sort_keys=True), encoding="utf-8")

    print(str(latest_md.resolve()))
    print(str(svg_path.resolve()))
    print(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
