"""Shared percolation analysis (Monkey B & Monkey C)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from scipy import stats
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

CONDITIONS = [
    ("Wakefulness", "MUA_W1", "LFP_W1"),
    ("Anaesthesia", "MUA_SO", "LFP_SO"),
    ("Awakening", "MUA_Aw", "LFP_Aw"),
]

PAIRWISE_LABELS = [
    ("Wakefulness", "Anaesthesia"),
    ("Anaesthesia", "Awakening"),
    ("Wakefulness", "Awakening"),
]


@dataclass
class SubjectConfig:
    name: str
    base: Path
    bad_ch: set[int]
    prune_threshold: float
    use_tanh: bool
    cross_sizes: list[int]
    mat_suffix: str = ".mat"
    save_figures: bool = True

    @property
    def ch(self) -> np.ndarray:
        return np.array([i for i in range(1, 193) if i not in self.bad_ch]) - 1

    @property
    def fig_dir(self) -> Path:
        d = self.base / "figures"
        d.mkdir(exist_ok=True)
        return d


SCRIPTS_ROOT = Path(__file__).resolve().parent.parent

MONKEY_B_CONFIG = SubjectConfig(
    name="Monkey B",
    base=SCRIPTS_ROOT / "MonkeyB" / "Data" / "Percolation",
    bad_ch={66, 20, 88, 174, 113, 131, 133},
    prune_threshold=0.5,
    use_tanh=True,
    cross_sizes=list(range(80, 125, 5)),
)

MONKEY_C_CONFIG = SubjectConfig(
    name="Monkey C",
    base=SCRIPTS_ROOT / "MonkeyC" / "Data" / "Percolation",
    bad_ch={62, 50},
    prune_threshold=0.1,
    use_tanh=False,
    cross_sizes=list(range(120, 165, 5)),
    mat_suffix="_C.mat",
)


def percolation(matrix: np.ndarray):
    """Match Percolation.m."""
    sort_value = np.unique(np.sort(matrix))
    matrix_corr = matrix.copy()
    matrix_work = matrix.copy()

    perc_threshold = None
    perc_matrix = None
    perc_threshold_step = None
    dont_enter = 0
    n_com_size = np.zeros((len(sort_value), 1))

    for t, val in enumerate(sort_value):
        matrix_work[matrix_work == val] = 0
        n_com, labels = connected_components(
            csr_matrix(np.triu(matrix_work, k=1)), directed=False
        )
        _, counts = np.unique(labels, return_counts=True)
        n_com_size[t, 0] = len(counts)

        if n_com > 1 and dont_enter == 0:
            dont_enter = 1
            matrix2 = matrix_corr.copy()
            corr_sort = sort_value[t - 1] if t > 0 else sort_value[0]
            matrix2[matrix2 <= corr_sort] = 0
            perc_threshold = corr_sort
            perc_matrix = matrix2
            perc_threshold_step = t - 1

    return perc_threshold, perc_matrix, sort_value, n_com_size, perc_threshold_step


def prune_matrix(perc_matrix: np.ndarray, threshold: float) -> np.ndarray:
    pruned = perc_matrix.copy()
    pruned[pruned < threshold] = 0
    return pruned


def analyze_condition(
    cfg: SubjectConfig, name: str, mua: np.ndarray, lfp: np.ndarray
):
    ind = round(mua.shape[1] / 16) - 1
    n_ch = len(cfg.ch)
    perc_stack = np.zeros((n_ch, n_ch, 16))
    sort_values = []
    n_com_sizes = []

    for i in range(16):
        ndx_start = ind * i
        ndx_end = ind + ind * i
        segment = mua[cfg.ch, ndx_start:ndx_end].T
        corr_mat = np.corrcoef(segment.T)
        _, perc_mat, sort_val, n_com_size, _ = percolation(corr_mat)
        if cfg.use_tanh:
            perc_stack[:, :, i] = np.tanh(perc_mat)
        else:
            perc_stack[:, :, i] = perc_mat
        sort_values.append(sort_val)
        n_com_sizes.append(n_com_size)

    if cfg.use_tanh:
        mean_perc = np.arctanh(np.mean(perc_stack, axis=2))
    else:
        mean_perc = np.mean(perc_stack, axis=2)

    pruned = prune_matrix(mean_perc, cfg.prune_threshold)

    return {
        "name": name,
        "ind": ind,
        "Perc_Matrix": mean_perc,
        "Perc_Matrix_pruned": pruned,
        "Sort_value": sort_values,
        "n_com_size": n_com_sizes,
    }


def cross_thresholds(result: dict, sizes: list[int]) -> np.ndarray:
    cross = np.full((16, len(sizes)), np.nan)
    for i in range(16):
        for n, j in enumerate(sizes):
            matches = np.where(result["n_com_size"][i].ravel() == j)[0]
            if len(matches):
                cross[i, n] = result["Sort_value"][i][matches[0]]
    return cross


def cross_threshold_per_size_std(cross: np.ndarray) -> np.ndarray:
    return np.nanstd(cross, axis=0)


def cross_threshold_statistics(
    cross_arrays: list[np.ndarray], labels: list[str]
) -> dict:
    groups = {
        label: cross_threshold_per_size_std(arr)
        for label, arr in zip(labels, cross_arrays)
    }
    arrays = [groups[label] for label in labels]
    if not all(len(a) == len(arrays[0]) for a in arrays):
        raise ValueError("Cross-threshold arrays must share the same component sizes.")

    stat, p = stats.friedmanchisquare(*arrays)
    result = {
        "n_component_sizes": len(arrays[0]),
        "omnibus": "Friedman",
        "omnibus_stat": float(stat),
        "omnibus_p": float(p),
        "per_size_std": groups,
        "condition_means": {
            label: float(np.mean(groups[label])) for label in labels
        },
        "pairwise": [],
    }

    n_pairs = len(PAIRWISE_LABELS)
    for a, b in PAIRWISE_LABELS:
        wstat, wp = stats.wilcoxon(groups[a], groups[b], alternative="two-sided")
        tstat, tp = stats.ttest_rel(groups[a], groups[b])
        result["pairwise"].append(
            {
                "comparison": f"{a} vs {b}",
                "wilcoxon_stat": float(wstat),
                "wilcoxon_p": float(wp),
                "wilcoxon_p_bonferroni": float(min(wp * n_pairs, 1.0)),
                "ttest_stat": float(tstat),
                "ttest_p": float(tp),
                "ttest_p_bonferroni": float(min(tp * n_pairs, 1.0)),
                "mean_diff": float(np.mean(groups[a] - groups[b])),
            }
        )
    return result


def _plot_heatmap(pruned: np.ndarray, title: str, out_path: Path):
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(pruned, cmap="viridis", aspect="auto")
    ax.set_title(title)
    plt.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _plot_graph(pruned: np.ndarray, title: str, out_path: Path):
    n = pruned.shape[0]
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pos = np.column_stack([np.cos(theta), np.sin(theta)])
    node_colors = np.array(["red" if i < 93 else "green" for i in range(n)])

    fig, ax = plt.subplots(figsize=(10, 10))
    rows, cols = np.where(np.triu(pruned, k=1) > 0)
    for r, c in zip(rows, cols):
        ax.plot(
            [pos[r, 0], pos[c, 0]],
            [pos[r, 1], pos[c, 1]],
            "--",
            color="gray",
            alpha=0.2,
            linewidth=0.5,
        )
    ax.scatter(pos[:, 0], pos[:, 1], c=node_colors, s=15, zorder=3)
    ax.set_title(title)
    ax.axis("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def run_subject(cfg: SubjectConfig) -> dict:
    """Run full percolation pipeline for one animal."""
    state_files = {
        "Wakefulness": f"Wakefulness{cfg.mat_suffix}",
        "Anaesthesia": f"Anaesthesia{cfg.mat_suffix}",
        "Awakening": f"Awakening{cfg.mat_suffix}",
    }

    conditions = {}
    for state, mat_name in state_files.items():
        mua_key, lfp_key = CONDITIONS[[c[0] for c in CONDITIONS].index(state)][1:]
        data = scipy.io.loadmat(cfg.base / mat_name)
        label = state
        print(f"\n=== {cfg.name} – {label} ===")
        print(
            f"  MUA shape: {data[mua_key].shape}, "
            f"segment length ind={round(data[mua_key].shape[1]/16)-1}"
        )
        conditions[label] = analyze_condition(cfg, label, data[mua_key], data[lfp_key])
        pruned = conditions[label]["Perc_Matrix_pruned"]
        n_edges = int(np.sum(pruned > 0) / 2)
        print(f"  Channels used: {pruned.shape[0]}")
        print(f"  Pruned matrix: {n_edges} edges, max weight {pruned.max():.4f}")
        print(
            f"  Fraction of pairs above {cfg.prune_threshold}: "
            f"{np.mean(pruned > 0):.4f}"
        )

        if cfg.save_figures:
            stem = mat_name.replace(".mat", "")
            _plot_heatmap(
                pruned,
                f"{cfg.name} {label} – pruned percolation matrix",
                cfg.fig_dir / f"percolation_heatmap_{stem}.png",
            )
            _plot_graph(
                pruned,
                f"{cfg.name} {label} – graph",
                cfg.fig_dir / f"percolation_graph_{stem}.png",
            )

    if cfg.save_figures:
        fig, ax = plt.subplots(figsize=(10, 6))
        colors = {"Wakefulness": "red", "Anaesthesia": "blue", "Awakening": "cyan"}
        for label, res in conditions.items():
            for i in range(16):
                ax.plot(
                    res["Sort_value"][i],
                    res["n_com_size"][i],
                    color=colors[label],
                    alpha=0.35,
                    linewidth=0.8,
                )
        ax.set_xlabel("Correlation threshold (sorted link values)")
        ax.set_ylabel("Number of connected components")
        ax.set_title(f"{cfg.name} – percolation curves (16 segments × 3 conditions)")
        from matplotlib.lines import Line2D

        ax.legend(
            handles=[
                Line2D([0], [0], color=c, lw=2, label=l) for l, c in colors.items()
            ],
            loc="upper right",
        )
        fig.tight_layout()
        fig.savefig(cfg.fig_dir / "percolation_curves_overlay.png", dpi=150)
        plt.close(fig)

    labels_short = [c[0] for c in CONDITIONS]
    cross = {
        label: cross_thresholds(res, cfg.cross_sizes)
        for label, res in conditions.items()
    }
    cross_arrays = [cross[label] for label in labels_short]

    x = [np.mean(np.nanstd(arr, axis=0)) for arr in cross_arrays]
    yerr = [np.std(np.nanstd(arr, axis=0)) for arr in cross_arrays]

    if cfg.save_figures:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.errorbar(range(1, 4), x, yerr=yerr, fmt="o-", capsize=5)
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(labels_short)
        ax.set_xlim(0.5, 3.5)
        ax.set_xlabel("Condition")
        size_label = f"{cfg.cross_sizes[0]}:{cfg.cross_sizes[1]-cfg.cross_sizes[0]}:{cfg.cross_sizes[-1]}"
        ax.set_ylabel(f"mean(std(threshold @ component size {size_label}))")
        ax.set_title(f"{cfg.name} – cross-threshold summary")
        fig.tight_layout()
        fig.savefig(cfg.fig_dir / "percolation_cross_threshold_errorbar.png", dpi=150)
        plt.close(fig)

    print(f"\n=== {cfg.name} cross-threshold summary (sizes {cfg.cross_sizes[0]}:"
          f"{cfg.cross_sizes[1]-cfg.cross_sizes[0]}:{cfg.cross_sizes[-1]}) ===")
    for lbl, arr in zip(labels_short, cross_arrays):
        print(
            f"  {lbl}: mean(std across segments)={np.nanmean(np.nanstd(arr, axis=0)):.4f}, "
            f"std(std)={np.std(np.nanstd(arr, axis=0)):.4f}"
        )

    stats_result = cross_threshold_statistics(cross_arrays, labels_short)
    print_cross_threshold_statistics(cfg.name, stats_result)

    if cfg.save_figures:
        print(f"  Figures saved to {cfg.fig_dir}")

    return {
        "config": cfg,
        "conditions": conditions,
        "cross": cross,
        "stats": stats_result,
        "errorbar_y": dict(zip(labels_short, x)),
    }


def print_cross_threshold_statistics(subject: str, stats_result: dict):
    n = stats_result["n_component_sizes"]
    print(f"\n=== {subject} cross-threshold statistics (n={n} component sizes) ===")
    print(
        f"  Friedman: chi2={stats_result['omnibus_stat']:.4f}, "
        f"p={stats_result['omnibus_p']:.4e}"
    )
    print("  Pairwise Wilcoxon (Bonferroni ×3):")
    for row in stats_result["pairwise"]:
        sig = "*" if row["wilcoxon_p_bonferroni"] < 0.05 else ""
        print(
            f"    {row['comparison']}: p_bonf={row['wilcoxon_p_bonferroni']:.4e}{sig}, "
            f"mean Δ={row['mean_diff']:.6f}, "
            f"t-test p_bonf={row['ttest_p_bonferroni']:.4e}"
        )


def summary_table_rows(subject: str, stats_result: dict) -> list[dict]:
    rows = []
    for row in stats_result["pairwise"]:
        rows.append(
            {
                "subject": subject,
                "comparison": row["comparison"],
                "mean_diff": row["mean_diff"],
                "wilcoxon_p": row["wilcoxon_p"],
                "wilcoxon_p_bonf": row["wilcoxon_p_bonferroni"],
                "ttest_p_bonf": row["ttest_p_bonferroni"],
            }
        )
    rows.insert(
        0,
        {
            "subject": subject,
            "comparison": "Friedman (omnibus)",
            "mean_diff": np.nan,
            "wilcoxon_p": stats_result["omnibus_p"],
            "wilcoxon_p_bonf": stats_result["omnibus_p"],
            "ttest_p_bonf": np.nan,
        },
    )
    for state, mean_val in stats_result["condition_means"].items():
        rows.append(
            {
                "subject": subject,
                "comparison": f"{state} (errorbar y)",
                "mean_diff": mean_val,
                "wilcoxon_p": np.nan,
                "wilcoxon_p_bonf": np.nan,
                "ttest_p_bonf": np.nan,
            }
        )
    return rows


def print_summary_table(all_rows: list[dict]):
    print("\n" + "=" * 100)
    print("SUMMARY TABLE – both animals")
    print("=" * 100)
    header = (
        f"{'Subject':<10} {'Comparison':<28} {'mean/std or Δ':>14} "
        f"{'Wilcoxon p':>12} {'p bonf':>12} {'t-test p bonf':>14}"
    )
    print(header)
    print("-" * len(header))
    for r in all_rows:
        md = r["mean_diff"]
        md_str = f"{md:.6f}" if not np.isnan(md) else ""
        wp = r["wilcoxon_p"]
        wp_str = f"{wp:.4e}" if not np.isnan(wp) else ""
        wpb = r["wilcoxon_p_bonf"]
        wpb_str = f"{wpb:.4e}" if not np.isnan(wpb) else ""
        tp = r["ttest_p_bonf"]
        tp_str = f"{tp:.4e}" if not np.isnan(tp) else ""
        print(
            f"{r['subject']:<10} {r['comparison']:<28} {md_str:>14} "
            f"{wp_str:>12} {wpb_str:>12} {tp_str:>14}"
        )
    print("=" * 100)
