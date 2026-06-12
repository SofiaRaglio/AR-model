#!/usr/bin/env python3
"""Python port of BoxplotInterArea.m — inter-area kernel performance across states."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Reuse loaders and tests from intra-area boxplot script
sys.path.insert(0, str(Path(__file__).resolve().parent))
import boxplot_kernels as bk  # noqa: E402

SCRIPTS_ROOT = bk.SCRIPTS_ROOT

# PM→PF (black in MATLAB), PF→PM (red)
DIRECTIONS = ("PM2PF", "PF2PM")
DIRECTION_COLORS = {"PM2PF": "black", "PF2PM": "#d62728"}


def kernel_path_interarea(
    subject: str, direction: str, transient: str, scripts_root: Path = SCRIPTS_ROOT
) -> Path:
    fname = f"KernelAndPerf{direction}_{transient}.mat"
    return scripts_root / subject / "Data" / "Kernel" / fname


def load_direction_states(
    subject: str, direction: str, scripts_root: Path = SCRIPTS_ROOT
) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for state, transient in bk.STATE_FILES[subject].items():
        path = kernel_path_interarea(subject, direction, transient, scripts_root)
        if not path.is_file():
            raise FileNotFoundError(path)
        out[state] = bk.load_performance(path)
    return out


def plot_interarea_subject(
    subject: str,
    data: dict[str, dict[str, np.ndarray]],
    out_path: Path,
) -> None:
    """One figure per monkey: PM2PF (black) and PF2PM (red), like BoxplotInterArea.m."""
    fig, ax = plt.subplots(figsize=(8, 5))
    n_states = len(bk.STATE_ORDER)
    width = 0.35

    for di, direction in enumerate(DIRECTIONS):
        offset = (di - 0.5) * width
        positions = np.arange(1, n_states + 1) + offset
        values = [data[direction][s] for s in bk.STATE_ORDER]
        bp = ax.boxplot(
            values,
            positions=positions,
            widths=width * 0.85,
            patch_artist=True,
            notch=True,
            showfliers=True,
        )
        color = DIRECTION_COLORS[direction]
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.55 if direction == "PM2PF" else 0.45)
            patch.set_edgecolor(color)

    ax.set_xticks(np.arange(1, n_states + 1))
    ax.set_xticklabels(bk.STATE_ORDER)
    ax.set_ylabel("PerformanceTest")
    label = bk.SUBJECT_LABELS.get(subject, subject)
    ax.set_title(f"{label} — inter-area kernel performance")
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=DIRECTION_COLORS[d], alpha=0.6)
        for d in DIRECTIONS
    ]
    ax.legend(handles, ["PM → PF", "PF → PM"], loc="best")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def format_p(p: float) -> str:
    return f"{p:.2e}" if p < 1e-3 else f"{p:.4g}"


def print_pairwise_table(rows: list[tuple]) -> None:
    print("\nPairwise p-values (Wilcoxon signed-rank, Bonferroni ×3, n=96 channels)")
    print("Monkey\tDirection\tWake vs Anest\tAnest vs Awak\tWake vs Awak")
    for monkey, direction, p_wa, p_aa, p_waw in rows:
        print(f"{monkey}\t{direction}\t{format_p(p_wa)}\t{format_p(p_aa)}\t{format_p(p_waw)}")


def run(scripts_root: Path, out_dir: Path, paired: bool = True) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    table_rows: list[tuple] = []

    for subject in bk.STATE_FILES:
        data = {
            direction: load_direction_states(subject, direction, scripts_root)
            for direction in DIRECTIONS
        }
        fig_path = out_dir / f"boxplot_interarea_{subject}.png"
        plot_interarea_subject(subject, data, fig_path)
        print(f"Saved: {fig_path}")

        for direction in DIRECTIONS:
            groups = data[direction]
            res = bk.state_significance(groups, paired=paired)
            label = bk.SUBJECT_LABELS.get(subject, subject)
            bk.print_significance(f"{label} — {direction}", res)

            pmap: dict[str, float] = {}
            for row in res["pairwise"]:
                if row["comparison"] == "Wakefulness vs Anesthesia":
                    pmap["wa_an"] = row["p_bonferroni"]
                elif row["comparison"] == "Anesthesia vs Awakening":
                    pmap["an_aw"] = row["p_bonferroni"]
                elif row["comparison"] == "Wakefulness vs Awakening":
                    pmap["wa_aw"] = row["p_bonferroni"]
            table_rows.append(
                (
                    bk.SUBJECT_LABELS.get(subject, subject),
                    direction,
                    pmap["wa_an"],
                    pmap["an_aw"],
                    pmap["wa_aw"],
                )
            )

    print_pairwise_table(table_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scripts-root",
        type=Path,
        default=SCRIPTS_ROOT,
        help="Scripts/ directory containing MonkeyB/ and MonkeyC/",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
    )
    parser.add_argument(
        "--unpaired",
        action="store_true",
        help="Use Kruskal-Wallis / Mann-Whitney instead of Friedman / Wilcoxon",
    )
    args = parser.parse_args()
    run(args.scripts_root, args.out_dir, paired=not args.unpaired)


if __name__ == "__main__":
    main()
