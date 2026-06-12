#!/usr/bin/env python3
"""Python port of Boxplot.m — kernel PerformanceTest across behavioral states."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from scipy import stats

SCRIPTS_ROOT = Path(__file__).resolve().parent.parent

STATE_FILES = {
    "MonkeyB": {
        "Wakefulness": "trans100",
        "Anesthesia": "trans3000",
        "Awakening": "trans7600",
    },
    "MonkeyC": {
        "Wakefulness": "trans1000",
        "Anesthesia": "trans4000",
        "Awakening": "trans6250",
    },
}

SUBJECT_LABELS = {
    "MonkeyB": "Monkey B",
    "MonkeyC": "Monkey C",
}

STATE_ORDER = ("Wakefulness", "Anesthesia", "Awakening")
PAIRWISE = (
    ("Wakefulness", "Anesthesia"),
    ("Wakefulness", "Awakening"),
    ("Anesthesia", "Awakening"),
)


def load_performance(mat_path: Path) -> np.ndarray:
    data = scipy.io.loadmat(mat_path)
    perf = np.asarray(data["PerformanceTest"]).squeeze()
    return perf.ravel()


def kernel_path(
    subject: str, area: str, transient: str, scripts_root: Path = SCRIPTS_ROOT
) -> Path:
    fname = f"KernelAndPerf{area}_{transient}.mat"
    return scripts_root / subject / "Data" / "Kernel" / fname


def load_subject_states(
    subject: str, area: str, scripts_root: Path = SCRIPTS_ROOT
) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for state, transient in STATE_FILES[subject].items():
        path = kernel_path(subject, area, transient, scripts_root)
        if not path.is_file():
            raise FileNotFoundError(path)
        out[state] = load_performance(path)
    return out


def state_significance(
    groups: dict[str, np.ndarray], paired: bool = True
) -> dict:
    """
    Test differences across Wakefulness / Anesthesia / Awakening.

    Default: Friedman (omnibus) + Wilcoxon signed-rank (pairwise, Bonferroni),
    treating the 96 channels as paired observations (same as channel-wise boxplots).
    """
    arrays = [groups[s] for s in STATE_ORDER]
    n = len(arrays[0])
    if not all(len(a) == n for a in arrays):
        raise ValueError("State vectors must have the same length (channels).")

    result: dict = {"n_channels": n, "paired": paired}

    if paired:
        stat, p = stats.friedmanchisquare(*arrays)
        result["omnibus"] = "Friedman"
        result["omnibus_stat"] = float(stat)
        result["omnibus_p"] = float(p)
        pairwise = []
        for a, b in PAIRWISE:
            wstat, wp = stats.wilcoxon(groups[a], groups[b], alternative="two-sided")
            pairwise.append(
                {
                    "comparison": f"{a} vs {b}",
                    "test": "Wilcoxon signed-rank",
                    "stat": float(wstat),
                    "p": float(wp),
                }
            )
        n_pairs = len(pairwise)
        for row in pairwise:
            row["p_bonferroni"] = min(row["p"] * n_pairs, 1.0)
        result["pairwise"] = pairwise
    else:
        stat, p = stats.kruskal(*arrays)
        result["omnibus"] = "Kruskal-Wallis"
        result["omnibus_stat"] = float(stat)
        result["omnibus_p"] = float(p)
        pairwise = []
        for a, b in PAIRWISE:
            ustat, up = stats.mannwhitneyu(groups[a], groups[b], alternative="two-sided")
            pairwise.append(
                {
                    "comparison": f"{a} vs {b}",
                    "test": "Mann-Whitney U",
                    "stat": float(ustat),
                    "p": float(up),
                }
            )
        n_pairs = len(pairwise)
        for row in pairwise:
            row["p_bonferroni"] = min(row["p"] * n_pairs, 1.0)
        result["pairwise"] = pairwise

    return result


def print_significance(label: str, res: dict) -> None:
    print(f"\n{'=' * 60}")
    print(label)
    print(f"  n channels = {res['n_channels']}, paired = {res['paired']}")
    print(
        f"  Omnibus ({res['omnibus']}): stat = {res['omnibus_stat']:.4f}, "
        f"p = {res['omnibus_p']:.4g}"
    )
    for row in res["pairwise"]:
        sig = "*" if row["p_bonferroni"] < 0.05 else "ns"
        print(
            f"  {row['comparison']}: {row['test']}, "
            f"p = {row['p']:.4g}, p_Bonferroni = {row['p_bonferroni']:.4g} {sig}"
        )


def plot_boxplots(
    data_by_subject: dict[str, dict[str, np.ndarray]],
    area: str,
    out_path: Path,
) -> None:
    """Two-subject overlay like Boxplot.m (notched boxplots)."""
    fig, ax = plt.subplots(figsize=(9, 5))
    subjects = list(data_by_subject.keys())
    n_states = len(STATE_ORDER)
    width = 0.35
    colors = {"MonkeyB": "#4daf4a", "MonkeyC": "#984ea3"}

    for si, subject in enumerate(subjects):
        offset = (si - 0.5) * width
        positions = np.arange(1, n_states + 1) + offset
        values = [data_by_subject[subject][s] for s in STATE_ORDER]
        bp = ax.boxplot(
            values,
            positions=positions,
            widths=width * 0.85,
            patch_artist=True,
            notch=True,
            showfliers=True,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(colors.get(subject, "gray"))
            patch.set_alpha(0.65)

    ax.set_xticks(np.arange(1, n_states + 1))
    ax.set_xticklabels(STATE_ORDER)
    ax.set_ylabel("PerformanceTest")
    ax.set_title(f"{area} cortex — kernel performance by state")
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=colors[s], alpha=0.65) for s in subjects
    ]
    ax.legend(
        handles, [SUBJECT_LABELS.get(s, s) for s in subjects], loc="best"
    )
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def run(scripts_root: Path, out_dir: Path, paired: bool) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    for area in ("PF", "PM"):
        data_by_subject: dict[str, dict[str, np.ndarray]] = {}
        for subject in STATE_FILES:
            data_by_subject[subject] = load_subject_states(
                subject, area, scripts_root
            )

        fig_path = out_dir / f"boxplot_{area}_performance.png"
        plot_boxplots(data_by_subject, area, fig_path)
        print(f"Saved: {fig_path}")

        for subject, groups in data_by_subject.items():
            res = state_significance(groups, paired=paired)
            label = SUBJECT_LABELS.get(subject, subject)
            print_significance(f"{label} — {area}", res)

        # Optional: omnibus across both animals (stack subjects — 192 "samples")
        stacked = {
            state: np.concatenate(
                [
                    data_by_subject["MonkeyB"][state],
                    data_by_subject["MonkeyC"][state],
                ]
            )
            for state in STATE_ORDER
        }
        res_pool = state_significance(stacked, paired=False)
        print_significance(f"Monkey B + Monkey C pooled (unpaired) — {area}", res_pool)


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
        help="Where to save figures",
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
