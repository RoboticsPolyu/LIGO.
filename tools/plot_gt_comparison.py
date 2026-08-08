#!/usr/bin/env python3
"""Plot and summarize the ECEF ground-truth comparison produced by LIGO."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
except ImportError as error:
    raise SystemExit(
        "Missing plotting dependency. Run: "
        "python3 -m pip install -r tools/requirements-plot.txt"
    ) from error


NUMERIC_COLUMNS = {
    "estimate_time", "gt_time", "match_dt", "gt_x", "gt_y", "gt_z",
    "est_x", "est_y", "est_z", "err_x", "err_y", "err_z", "err_3d",
    "err_horizontal", "err_vertical", "gt_vx", "gt_vy", "gt_vz",
    "est_vx", "est_vy", "est_vz", "vel_err_3d", "gt_quality",
    "satellites", "frame", "lambda_ratio",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot LIGO ECEF trajectory and ground-truth errors.")
    parser.add_argument(
        "input", nargs="?", default="Log/gt_comparison_ecef.txt",
        help="comparison text file (default: Log/gt_comparison_ecef.txt)")
    parser.add_argument("--output-dir", type=Path,
                        help="save plots and summary here")
    parser.add_argument("--start", type=float, default=0.0,
                        help="seconds after the first estimate to include")
    parser.add_argument("--end", type=float,
                        help="seconds after the first estimate to stop")
    parser.add_argument("--max-match-dt", type=float,
                        help="discard pairs with larger absolute time mismatch")
    parser.add_argument("--max-error", type=float,
                        help="discard samples whose 3D position error exceeds this [m]")
    parser.add_argument("--percentile", type=float, default=95.0,
                        help="percentile reported in statistics (default: 95)")
    parser.add_argument("--title", default="LIGO ground-truth evaluation")
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument("--time-cmap", default="viridis_r",
                        help="Matplotlib colormap for trajectory time (default: viridis_r)")
    parser.add_argument("--min-alpha", type=float, default=0.2,
                        help="opacity at trajectory start, in [0, 1] (default: 0.2)")
    parser.add_argument("--no-show", action="store_true",
                        help="do not open interactive plot windows")
    return parser.parse_args()


def load_report(path: Path) -> Dict[str, np.ndarray]:
    header: List[str] | None = None
    rows: List[List[str]] = []
    with path.open("r", encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                header = line[1:].split()
            elif header is not None:
                fields = line.split()
                if len(fields) == len(header):
                    rows.append(fields)
    if header is None:
        raise ValueError("missing '# ...' column header")
    if not rows:
        raise ValueError("report contains no matched samples")

    result: Dict[str, np.ndarray] = {}
    columns = list(zip(*rows))
    for name, values in zip(header, columns):
        if name in NUMERIC_COLUMNS:
            result[name] = np.asarray(values, dtype=float)
        else:
            result[name] = np.asarray(values, dtype=str)
    return result


def apply_filters(data: Dict[str, np.ndarray], args: argparse.Namespace) -> Dict[str, np.ndarray]:
    elapsed = data["estimate_time"] - data["estimate_time"][0]
    mask = np.isfinite(data["err_3d"]) & (elapsed >= args.start)
    if args.end is not None:
        mask &= elapsed <= args.end
    if args.max_match_dt is not None:
        mask &= np.abs(data["match_dt"]) <= args.max_match_dt
    if args.max_error is not None:
        mask &= data["err_3d"] <= args.max_error
    if not np.any(mask):
        raise ValueError("all samples were removed by the selected filters")
    filtered = {name: values[mask] for name, values in data.items()}
    filtered["elapsed"] = elapsed[mask]
    return filtered


def finite_stats(values: np.ndarray, percentile: float) -> tuple[float, float, float, float]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return (math.nan,) * 4
    return (float(np.sqrt(np.mean(values * values))), float(np.mean(values)),
            float(np.median(values)), float(np.percentile(values, percentile)))


def summary_text(data: Dict[str, np.ndarray], percentile: float) -> str:
    lines = [
        f"samples: {len(data['estimate_time'])}",
        f"duration_s: {data['elapsed'][-1] - data['elapsed'][0]:.3f}",
        f"mean_abs_match_dt_s: {np.mean(np.abs(data['match_dt'])):.6f}",
        "metric rmse mean median p{:.1f}".format(percentile),
    ]
    metrics = {
        "error_3d_m": data["err_3d"],
        "horizontal_m": data["err_horizontal"],
        "vertical_abs_m": np.abs(data["err_vertical"]),
        "velocity_3d_mps": data["vel_err_3d"],
    }
    for name, values in metrics.items():
        stats = finite_stats(values, percentile)
        lines.append(name + " " + " ".join(f"{value:.6f}" for value in stats))
    if "rtk_status" in data:
        names, counts = np.unique(data["rtk_status"], return_counts=True)
        lines.append("rtk_status: " + ", ".join(
            f"{name}={count}" for name, count in zip(names, counts)))
    return "\n".join(lines) + "\n"


def temporal_line_3d(ax, points: np.ndarray, elapsed: np.ndarray,
                     cmap: str, min_alpha: float, linewidth: float,
                     label: str) -> mpl.cm.ScalarMappable:
    """Add a 3D line whose color and opacity encode segment time."""
    segments = np.stack((points[:-1], points[1:]), axis=1)
    segment_time = 0.5 * (elapsed[:-1] + elapsed[1:])
    time_min, time_max = float(elapsed[0]), float(elapsed[-1])
    if time_max <= time_min:
        time_max = time_min + 1.0
    normalization = mpl.colors.Normalize(vmin=time_min, vmax=time_max)
    color_map = mpl.colormaps[cmap]
    colors = color_map(normalization(segment_time))
    progress = normalization(segment_time)
    colors[:, 3] = min_alpha + (1.0 - min_alpha) * progress
    collection = Line3DCollection(
        segments, colors=colors, linewidths=linewidth, label=label)
    ax.add_collection3d(collection)
    ax.auto_scale_xyz(points[:, 0], points[:, 1], points[:, 2])
    return mpl.cm.ScalarMappable(norm=normalization, cmap=color_map)


def make_plots(data: Dict[str, np.ndarray], title: str,
               time_cmap: str, min_alpha: float):
    t = data["elapsed"]
    gt = np.column_stack((data["gt_x"], data["gt_y"], data["gt_z"]))
    est = np.column_stack((data["est_x"], data["est_y"], data["est_z"]))
    origin = gt[0]

    trajectory = plt.figure(figsize=(13, 6), constrained_layout=True)
    trajectory.suptitle(title + " — trajectory")
    ax_xy = trajectory.add_subplot(1, 2, 1)
    ax_xy.plot(gt[:, 0] - origin[0], gt[:, 1] - origin[1], label="GT", linewidth=2)
    ax_xy.plot(est[:, 0] - origin[0], est[:, 1] - origin[1], label="Estimate", alpha=.85)
    ax_xy.set(xlabel="ECEF X relative [m]", ylabel="ECEF Y relative [m]")
    ax_xy.axis("equal")
    ax_xy.grid(True, alpha=.3)
    ax_xy.legend()
    ax_3d = trajectory.add_subplot(1, 2, 2, projection="3d")
    gt_relative = gt - origin
    est_relative = est - origin
    time_map = temporal_line_3d(
        ax_3d, gt_relative, t, time_cmap, min_alpha, 3.0, "GT")
    temporal_line_3d(
        ax_3d, est_relative, t, time_cmap, min_alpha, 1.2, "Estimate")
    # Start/end markers make the direction explicit even in a printed figure.
    ax_3d.scatter(*gt_relative[0], color=time_map.cmap(0.0), marker="o", s=35)
    ax_3d.scatter(*gt_relative[-1], color=time_map.cmap(1.0), marker="s", s=35)
    ax_3d.set(xlabel="dX [m]", ylabel="dY [m]", zlabel="dZ [m]")
    ax_3d.legend()
    colorbar = trajectory.colorbar(time_map, ax=ax_3d, pad=0.12, shrink=0.75)
    colorbar.set_label("Elapsed time [s] (lighter → earlier, darker → later)")

    errors, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    errors.suptitle(title + " — errors")
    axes[0, 0].plot(t, data["err_x"], label="X")
    axes[0, 0].plot(t, data["err_y"], label="Y")
    axes[0, 0].plot(t, data["err_z"], label="Z")
    axes[0, 0].set(ylabel="ECEF error [m]")
    axes[0, 0].legend(ncol=3)
    axes[0, 1].plot(t, data["err_3d"], label="3D")
    axes[0, 1].plot(t, data["err_horizontal"], label="Horizontal")
    axes[0, 1].plot(t, np.abs(data["err_vertical"]), label="|Vertical|")
    axes[0, 1].set(ylabel="Position error [m]")
    axes[0, 1].legend()
    axes[1, 0].plot(t, data["vel_err_3d"], color="tab:purple")
    axes[1, 0].set(xlabel="Elapsed time [s]", ylabel="3D velocity error [m/s]")
    axes[1, 1].plot(t, np.abs(data["match_dt"]) * 1000, label="Time match")
    ax_sat = axes[1, 1].twinx()
    ax_sat.plot(t, data["satellites"], color="tab:green", alpha=.55, label="Satellites")
    axes[1, 1].set(xlabel="Elapsed time [s]", ylabel="|match dt| [ms]")
    ax_sat.set_ylabel("Satellites")
    for axis in axes.flat:
        axis.grid(True, alpha=.3)
    return trajectory, errors


def main() -> None:
    args = parse_args()
    path = Path(args.input).expanduser()
    data = apply_filters(load_report(path), args)
    text = summary_text(data, args.percentile)
    print(text, end="")
    if not 0.0 <= args.min_alpha <= 1.0:
        raise ValueError("--min-alpha must be between 0 and 1")
    figures = make_plots(data, args.title, args.time_cmap, args.min_alpha)
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        figures[0].savefig(args.output_dir / "trajectory_ecef.png", dpi=args.dpi)
        figures[1].savefig(args.output_dir / "errors.png", dpi=args.dpi)
        (args.output_dir / "summary.txt").write_text(text, encoding="utf-8")
        print(f"saved results to {args.output_dir}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
