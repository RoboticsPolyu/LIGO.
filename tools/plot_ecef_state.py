#!/usr/bin/env python3
"""Plot the estimate-only ECEF state log produced by LIGO."""

from __future__ import annotations

import argparse
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


STRING_COLUMNS = {"rtk_status"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot LIGO ECEF state output.")
    parser.add_argument(
        "input", nargs="?", default="Log/state_estimate_ecef.txt",
        help="ECEF state text file (default: Log/state_estimate_ecef.txt)")
    parser.add_argument("--output-dir", type=Path,
                        help="save figures and summary to this directory")
    parser.add_argument("--start", type=float, default=0.0,
                        help="seconds after the first state to include")
    parser.add_argument("--end", type=float,
                        help="seconds after the first state to stop")
    parser.add_argument("--position", choices=("antenna", "imu"),
                        default="antenna",
                        help="position used for trajectory plots (default: antenna)")
    parser.add_argument("--time-cmap", default="viridis_r")
    parser.add_argument("--min-alpha", type=float, default=0.2)
    parser.add_argument("--title", default="LIGO ECEF state estimate")
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument("--no-show", action="store_true")
    return parser.parse_args()


def load_log(path: Path) -> Dict[str, np.ndarray]:
    header: List[str] | None = None
    rows: List[List[str]] = []
    with path.open("r", encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                header = line[1:].split()
                continue
            if header is not None:
                fields = line.split()
                if len(fields) == len(header):
                    rows.append(fields)
    if header is None:
        raise ValueError("missing '# ...' column header")
    if not rows:
        raise ValueError("ECEF state file contains no samples")
    columns = list(zip(*rows))
    return {
        name: np.asarray(values, dtype=str if name in STRING_COLUMNS else float)
        for name, values in zip(header, columns)
    }


def filter_time(data: Dict[str, np.ndarray], start: float,
                end: float | None) -> Dict[str, np.ndarray]:
    elapsed = data["time"] - data["time"][0]
    mask = np.isfinite(data["time"]) & (elapsed >= start)
    if end is not None:
        mask &= elapsed <= end
    if not np.any(mask):
        raise ValueError("selected time range contains no samples")
    result = {name: values[mask] for name, values in data.items()}
    result["elapsed"] = elapsed[mask]
    return result


def colored_3d_line(ax, points: np.ndarray, elapsed: np.ndarray,
                    cmap_name: str, min_alpha: float):
    if len(points) < 2:
        ax.scatter(*points[0], marker="o")
        return mpl.cm.ScalarMappable(
            norm=mpl.colors.Normalize(0.0, 1.0), cmap=mpl.colormaps[cmap_name])
    segments = np.stack((points[:-1], points[1:]), axis=1)
    middle_time = 0.5 * (elapsed[:-1] + elapsed[1:])
    time_min, time_max = float(elapsed[0]), float(elapsed[-1])
    if time_max <= time_min:
        time_max = time_min + 1.0
    normalization = mpl.colors.Normalize(time_min, time_max)
    color_map = mpl.colormaps[cmap_name]
    colors = color_map(normalization(middle_time))
    colors[:, 3] = min_alpha + (1.0 - min_alpha) * normalization(middle_time)
    ax.add_collection3d(Line3DCollection(segments, colors=colors, linewidths=2.0))
    ax.auto_scale_xyz(points[:, 0], points[:, 1], points[:, 2])
    return mpl.cm.ScalarMappable(norm=normalization, cmap=color_map)


def make_figures(data: Dict[str, np.ndarray], args: argparse.Namespace):
    t = data["elapsed"]
    prefix = "antenna" if args.position == "antenna" else "imu"
    position = np.column_stack(tuple(data[f"{prefix}_{axis}"] for axis in "xyz"))
    relative = position - position[0]

    trajectory = plt.figure(figsize=(13, 6), constrained_layout=True)
    trajectory.suptitle(args.title + f" — {args.position} trajectory")
    ax_xy = trajectory.add_subplot(1, 2, 1)
    scatter = ax_xy.scatter(relative[:, 0], relative[:, 1], c=t,
                            cmap=args.time_cmap, s=6, alpha=.8)
    ax_xy.plot(relative[:, 0], relative[:, 1], color="0.65", linewidth=.5)
    ax_xy.scatter(*relative[0, :2], marker="o", color=scatter.cmap(0.0), label="Start")
    ax_xy.scatter(*relative[-1, :2], marker="s", color=scatter.cmap(1.0), label="End")
    ax_xy.set(xlabel="relative ECEF X [m]", ylabel="relative ECEF Y [m]")
    ax_xy.axis("equal")
    ax_xy.grid(True, alpha=.3)
    ax_xy.legend()

    ax_3d = trajectory.add_subplot(1, 2, 2, projection="3d")
    time_map = colored_3d_line(ax_3d, relative, t, args.time_cmap, args.min_alpha)
    ax_3d.scatter(*relative[0], marker="o", color=time_map.cmap(0.0), s=35)
    ax_3d.scatter(*relative[-1], marker="s", color=time_map.cmap(1.0), s=35)
    ax_3d.set(xlabel="dX [m]", ylabel="dY [m]", zlabel="dZ [m]")
    colorbar = trajectory.colorbar(time_map, ax=ax_3d, pad=.12, shrink=.75)
    colorbar.set_label("Elapsed time [s]")

    state, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    state.suptitle(args.title + " — state history")
    for axis_name, color in zip("xyz", ("tab:red", "tab:green", "tab:blue")):
        axes[0, 0].plot(t, data[f"{prefix}_{axis_name}"] - data[f"{prefix}_{axis_name}"][0],
                        label=axis_name.upper(), color=color)
        axes[0, 1].plot(t, data[f"vel_{axis_name}"],
                        label=axis_name.upper(), color=color)
        axes[1, 0].plot(t, data[f"ba_{axis_name}"], linestyle="-",
                        label=f"ba {axis_name.upper()}", color=color)
        axes[1, 0].plot(t, data[f"bg_{axis_name}"], linestyle="--",
                        label=f"bg {axis_name.upper()}", color=color, alpha=.75)
    axes[0, 0].set(ylabel="Relative ECEF position [m]")
    axes[0, 1].set(ylabel="ECEF velocity [m/s]")
    axes[1, 0].set(xlabel="Elapsed time [s]", ylabel="Bias")
    axes[0, 0].legend(ncol=3)
    axes[0, 1].legend(ncol=3)
    axes[1, 0].legend(ncol=3, fontsize=8)

    axes[1, 1].plot(t, data["satellites"], color="tab:green", label="Satellites")
    axes[1, 1].set(xlabel="Elapsed time [s]", ylabel="Satellites")
    ratio_axis = axes[1, 1].twinx()
    ratio_axis.plot(t, data["lambda_ratio"], color="tab:purple", alpha=.7,
                    label="LAMBDA ratio")
    ratio_axis.set_ylabel("LAMBDA ratio")
    for axis in axes.flat:
        axis.grid(True, alpha=.3)
    return trajectory, state


def make_summary(data: Dict[str, np.ndarray]) -> str:
    speed = np.sqrt(data["vel_x"] ** 2 + data["vel_y"] ** 2 + data["vel_z"] ** 2)
    statuses, counts = np.unique(data["rtk_status"], return_counts=True)
    return "\n".join((
        f"samples: {len(data['time'])}",
        f"duration_s: {data['elapsed'][-1] - data['elapsed'][0]:.3f}",
        f"mean_speed_mps: {np.mean(speed):.6f}",
        f"max_speed_mps: {np.max(speed):.6f}",
        f"mean_satellites: {np.mean(data['satellites']):.3f}",
        "rtk_status: " + ", ".join(
            f"{name}={count}" for name, count in zip(statuses, counts)),
    )) + "\n"


def main() -> None:
    args = parse_args()
    if not 0.0 <= args.min_alpha <= 1.0:
        raise ValueError("--min-alpha must be between 0 and 1")
    data = filter_time(load_log(Path(args.input).expanduser()), args.start, args.end)
    summary = make_summary(data)
    print(summary, end="")
    figures = make_figures(data, args)
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        figures[0].savefig(args.output_dir / "ecef_trajectory.png", dpi=args.dpi)
        figures[1].savefig(args.output_dir / "ecef_state_history.png", dpi=args.dpi)
        (args.output_dir / "summary.txt").write_text(summary, encoding="utf-8")
        print(f"saved results to {args.output_dir}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
