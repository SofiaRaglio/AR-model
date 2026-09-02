#!/usr/bin/env python3
"""
ARMA LFP prediction: per-channel (channel-to-channel) and area-average.

Mirrors ARMA_PF.m / ARMA_PM.m for both animals and three behavioral stages
(Wakefulness, Anesthesia, Awakening) with subject-specific transient periods.

Comparisons (to justify per-channel analysis):
  1. Mean ARMA: fit/predict mean(LFP) from lagged mean(MUA).
  2. Mean of per-channel rho: average test correlation across channel-wise fits.
  3. Mean from per-channel preds: corr(mean_t y_hat_ch(t), mean_t y_ch(t)) using
     stored kernels — pooling reconstructions after per-channel fitting.
  4. Per-channel ARMA on LFP with the area mean subtracted (demeaned LFP).

Plots spatial maps, boxplots, and a direct three-way comparison of mean vs
per-channel summaries.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from scipy import stats
from scipy.linalg import qr, solve_triangular
from sklearn.decomposition import PCA

from evaluate_kernel_metrics import (
    KERNEL_LENGTH,
    SMOOTHING_WINDOW,
    area_config,
    kernel_to_alpha,
    load_kernel_mat,
    load_subject_data,
    select_mua_lfp,
    test_indices,
)

SCRIPTS_ROOT = Path(__file__).resolve().parent.parent

MONKEY_LABELS = {
    "MonkeyB": "Monkey B",
    "MonkeyC": "Monkey C",
}

# MEAMap layouts and channel lists from ARMA_PF.m / ARMA_PM.m
AREA_CONFIG: dict[str, dict[str, dict]] = {
    "MonkeyB": {
        "PF": {
            "ch": sorted(set(range(1, 97)) - {66, 20, 88}),
            "map": np.array(
                [
                    [0, 9, 8, 7, 10, 11, 12, 13, 15, 0],
                    [72, 73, 40, 41, 5, 4, 3, 14, 16, 17],
                    [71, 74, 39, 42, 6, 32, 2, 29, 18, 19],
                    [70, 75, 38, 43, 48, 1, 31, 28, 27, 20],
                    [69, 76, 37, 44, 45, 49, 51, 30, 26, 21],
                    [68, 77, 36, 35, 46, 47, 50, 55, 25, 22],
                    [67, 78, 34, 33, 63, 52, 53, 54, 56, 23],
                    [66, 79, 81, 64, 62, 61, 60, 59, 58, 24],
                    [65, 80, 82, 83, 93, 92, 91, 87, 57, 89],
                    [0, 96, 95, 94, 84, 85, 86, 90, 88, 0],
                ]
            ),
        },
        "PM": {
            "ch": sorted(set(range(97, 193)) - {174, 113, 131, 133}),
            "map": np.array(
                [
                    [0, 105, 104, 103, 106, 107, 108, 109, 0, 111],
                    [168, 169, 136, 137, 101, 100, 99, 110, 112, 113],
                    [167, 170, 135, 138, 102, 128, 98, 125, 114, 115],
                    [166, 171, 134, 139, 144, 97, 127, 124, 123, 116],
                    [165, 172, 133, 140, 141, 145, 147, 126, 122, 117],
                    [164, 173, 132, 131, 142, 143, 146, 151, 121, 118],
                    [163, 174, 130, 129, 159, 148, 149, 150, 152, 119],
                    [162, 175, 177, 160, 158, 157, 156, 155, 154, 120],
                    [161, 176, 178, 179, 189, 188, 187, 183, 153, 185],
                    [0, 192, 191, 190, 180, 181, 182, 186, 184, 0],
                ]
            ),
        },
    },
    "MonkeyC": {
        "PF": {
            "ch": sorted(set(range(1, 97)) - {50, 62}),
            "map": np.array(
                [
                    [0, 17, 19, 20, 21, 22, 23, 24, 89, 0],
                    [15, 16, 18, 27, 26, 25, 56, 58, 57, 88],
                    [13, 14, 29, 28, 30, 55, 54, 59, 87, 90],
                    [12, 3, 2, 31, 51, 50, 53, 60, 91, 86],
                    [11, 4, 32, 1, 49, 47, 52, 61, 92, 85],
                    [10, 5, 6, 48, 45, 46, 63, 62, 93, 84],
                    [7, 41, 42, 43, 44, 35, 33, 64, 83, 94],
                    [8, 40, 39, 38, 37, 36, 34, 81, 82, 95],
                    [9, 73, 74, 75, 76, 77, 78, 79, 80, 96],
                    [0, 72, 71, 70, 69, 68, 67, 66, 65, 0],
                ]
            ),
        },
        "PM": {
            "ch": list(range(97, 193)),
            "map": np.array(
                [
                    [0, 113, 115, 116, 117, 118, 119, 120, 185, 0],
                    [111, 112, 114, 123, 122, 121, 152, 154, 153, 184],
                    [109, 110, 125, 124, 126, 151, 150, 155, 183, 186],
                    [108, 99, 98, 127, 147, 146, 149, 156, 187, 182],
                    [107, 100, 128, 97, 145, 143, 148, 157, 188, 181],
                    [106, 101, 102, 144, 141, 142, 159, 158, 189, 180],
                    [103, 137, 138, 139, 140, 131, 129, 160, 179, 190],
                    [104, 136, 135, 134, 133, 132, 130, 177, 178, 191],
                    [105, 169, 170, 171, 172, 173, 174, 175, 176, 192],
                    [0, 168, 167, 166, 165, 164, 163, 162, 161, 0],
                ]
            ),
        },
    },
}

STATE_FILES = {
    "MonkeyB": {
        "Wakefulness": 100,
        "Anesthesia": 3000,
        "Awakening": 7600,
    },
    "MonkeyC": {
        "Wakefulness": 1000,
        "Anesthesia": 4000,
        "Awakening": 6250,
    },
}

STATE_ORDER = ("Wakefulness", "Anesthesia", "Awakening")
LAG_DT_S = 0.005
ANESTHESIA_STATE = "Anesthesia"


def kernel_path(subject: str, area: str, transient: int) -> Path:
    return SCRIPTS_ROOT / subject / "Kernel" / f"KernelAndPerf{area}_trans{transient}.mat"


def load_performance_test(mat_path: Path) -> np.ndarray:
    data = scipy.io.loadmat(mat_path, variable_names=["PerformanceTest"])
    return np.asarray(data["PerformanceTest"]).squeeze().astype(float)


def build_lagged_scalar(mua_1d: np.ndarray, sample_idx: np.ndarray) -> np.ndarray:
    """Lagged input from a single mean-MUA trace: (kernel_length, n_samples)."""
    offsets = np.arange(1, KERNEL_LENGTH + 1, SMOOTHING_WINDOW)
    return np.vstack([mua_1d[sample_idx - off] for off in offsets])


def build_lagged_input_fast(
    mua: np.ndarray, row_idx: np.ndarray, sample_idx: np.ndarray
) -> np.ndarray:
    """Stack lagged MUA blocks: (n_ch * kernel_length, n_samples)."""
    mua_sub = np.ascontiguousarray(mua[row_idx])
    offsets = np.arange(1, KERNEL_LENGTH + 1, SMOOTHING_WINDOW)
    n_ch = mua_sub.shape[0]
    n_samples = len(sample_idx)
    out = np.empty((len(offsets) * n_ch, n_samples), dtype=np.float64)
    for i, off in enumerate(offsets):
        out[i * n_ch : (i + 1) * n_ch] = mua_sub[:, sample_idx - off]
    return out


def arma_fit(
    y_train: np.ndarray,
    x_train: np.ndarray,
    y_test: np.ndarray,
    x_test: np.ndarray,
) -> tuple[float, float, float, float]:
    """Train linear ARMA readout; return train/test rho and p-values."""
    gram = x_train @ x_train.T
    alpha = np.linalg.solve(gram, x_train @ y_train)
    y_hat_train = alpha @ x_train
    y_hat_test = alpha @ x_test
    rho_tr, p_tr = stats.pearsonr(y_train, y_hat_train)
    rho_te, p_te = stats.pearsonr(y_test, y_hat_test)
    return float(rho_tr), float(p_tr), float(rho_te), float(p_te)


def fit_mean_kernel_alpha(
    mua: np.ndarray,
    lfp: np.ndarray,
    cfg,
    transient: float,
    time: np.ndarray,
) -> np.ndarray:
    """Fitted temporal kernel for mean LFP ← mean MUA (length KERNEL_LENGTH)."""
    ndx_lp, _ = test_indices(time, transient)
    mua_avg = mua[cfg.mua_rows].mean(axis=0)
    lfp_avg = lfp[cfg.lfp_rows].mean(axis=0)
    in_train = build_lagged_scalar(mua_avg, ndx_lp)
    y_tr = lfp_avg[ndx_lp + 1]
    gram = in_train @ in_train.T
    return np.linalg.solve(gram, in_train @ y_tr)


def average_arma(
    mua: np.ndarray,
    lfp: np.ndarray,
    cfg,
    transient: float,
    time: np.ndarray,
) -> tuple[float, float, float, float]:
    """Predict mean LFP from mean MUA (with temporal kernel)."""
    ndx_lp, ndx_tp = test_indices(time, transient)
    mua_avg = mua[cfg.mua_rows].mean(axis=0)
    lfp_avg = lfp[cfg.lfp_rows].mean(axis=0)
    in_train = build_lagged_scalar(mua_avg, ndx_lp)
    in_test = build_lagged_scalar(mua_avg, ndx_tp)
    y_tr = lfp_avg[ndx_lp + 1]
    y_te = lfp_avg[ndx_tp + 1]
    return arma_fit(y_tr, in_train, y_te, in_test)


def kernel_pca_pc1(kernel: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
    """PCA on kernels (MATLAB: Z = reshape(Kernel, NoT, NoC*NoCT); pca(Z'))."""
    no_t, no_c, no_ct = kernel.shape
    z = kernel.reshape(no_t, no_c * no_ct, order="F")
    pca = PCA(n_components=1)
    scores = pca.fit_transform(z.T)[:, 0]
    return pca.components_[0], float(pca.explained_variance_ratio_[0]), scores


def pc1_scores_per_output(kernel: np.ndarray, coeff_pc1: np.ndarray) -> np.ndarray:
    """PC1 projection norm per output channel (KernelDistDependence-style)."""
    no_ct = kernel.shape[2]
    out_scores = np.zeros(no_ct)
    for k in range(no_ct):
        proj = coeff_pc1 @ kernel[:, :, k]
        out_scores[k] = float(np.sqrt(np.mean(proj**2)))
    return out_scores


def demean_lfp_by_area(lfp: np.ndarray, cfg) -> np.ndarray:
    """Subtract area-mean LFP from each channel (common component removed)."""
    lfp_avg = lfp[cfg.lfp_rows].mean(axis=0)
    out = lfp.copy()
    for row in cfg.lfp_rows:
        out[row] = lfp[row] - lfp_avg
    return out


def mean_from_per_channel_predictions(
    kernels: np.ndarray,
    lfp: np.ndarray,
    cfg,
    in_test: np.ndarray,
    ndx_tp: np.ndarray,
) -> tuple[float, float]:
    """Correlation of area-mean LFP vs mean of per-channel kernel reconstructions."""
    y_true = lfp[cfg.lfp_rows][:, ndx_tp + 1].mean(axis=0)
    preds = []
    for k_idx in cfg.kernel_targets:
        alpha = kernel_to_alpha(kernels[:, :, k_idx])
        preds.append(alpha @ in_test)
    y_pred = np.mean(preds, axis=0)
    rho, p = stats.pearsonr(y_true, y_pred)
    return float(rho), float(p)


def per_channel_arma(
    lfp: np.ndarray,
    cfg,
    in_train: np.ndarray,
    in_test: np.ndarray,
    ndx_lp: np.ndarray,
    ndx_tp: np.ndarray,
    *,
    demean_lfp: bool = False,
) -> np.ndarray:
    """Fit ARMA per output channel; return PerformanceTest-style vector."""
    lfp_use = demean_lfp_by_area(lfp, cfg) if demean_lfp else lfp

    y_tr = lfp_use[cfg.lfp_rows][:, ndx_lp + 1]
    y_te = lfp_use[cfg.lfp_rows][:, ndx_tp + 1]
    q, r = qr(in_train.T, mode="economic")
    alphas = solve_triangular(r, q.T @ y_tr.T)
    y_hat_te = alphas.T @ in_test

    y_c = y_te - y_te.mean(axis=1, keepdims=True)
    p_c = y_hat_te - y_hat_te.mean(axis=1, keepdims=True)
    rho = (y_c * p_c).sum(axis=1) / np.sqrt(
        (y_c**2).sum(axis=1) * (p_c**2).sum(axis=1)
    )
    perf_test = np.full(lfp.shape[0], np.nan)
    for i, lfp_row in enumerate(cfg.lfp_rows):
        perf_test[lfp_row] = float(rho[i])
    return perf_test


def pf_all_ch() -> np.ndarray:
    return np.arange(1, 97)


def pm_all_ch() -> np.ndarray:
    return np.arange(97, 193)


def lfp_row_for_channel(ch: int, area: str) -> int:
    return ch - 1 if area == "PF" else ch - 97


def channel_values_from_perf(perf_test: np.ndarray, cfg) -> np.ndarray:
    return np.array([perf_test[r] for r in cfg.lfp_rows if np.isfinite(perf_test[r])])


def plot_state_row(
    results: dict[str, dict],
    mea_map: np.ndarray,
    ch_list: list[int],
    all_ch: np.ndarray,
    monkey: str,
    area: str,
    output: Path | None,
) -> plt.Figure:
    """Three-panel spatial map (one per behavioral state)."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    area_name = "Prefrontal" if area == "PF" else "Premotor"
    im = None

    for ax, state in zip(axes, STATE_ORDER):
        res = results[state]
        measure = np.zeros_like(mea_map, dtype=float)
        perf = res["perf_test"]
        for ch in all_ch:
            r, c = np.where(mea_map == ch)
            if r.size:
                row = lfp_row_for_channel(ch, area)
                measure[r[0], c[0]] = perf[row]
        im = ax.imshow(measure, vmin=0, vmax=0.8, cmap="RdYlBu_r", origin="upper")
        for ch in ch_list:
            r, c = np.where(mea_map == ch)
            if r.size:
                ax.text(c[0], r[0], str(ch), ha="center", va="center", fontsize=6)
        ax.set_title(
            f"{state}\n"
            f"per-ch median r = {np.median(res['perf_ch_values']):.3f}  |  "
            f"mean ARMA r = {res['rho_avg_test']:.3f}  |  "
            f"mean(pred) r = {res['rho_mean_pred']:.3f}",
            fontsize=8,
        )
        ax.set_xlabel("X MEA")
        ax.set_ylabel("Y MEA")

    if im is not None:
        cbar = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02)
        cbar.set_label("Test Pearson r (per channel)")
    fig.suptitle(
        f"{MONKEY_LABELS[monkey]} — {area_name}  |  "
        "per-channel ARMA (map) vs mean ARMA / mean(per-ch preds)",
        fontsize=12,
        y=1.03,
    )
    fig.tight_layout()
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=200, bbox_inches="tight")
        fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
        print(f"Saved {output}")
    return fig


MONKEY_METRIC_COLORS = {
    "MonkeyB": ("#1b7837", "#4daf4a", "#a6d96a"),
    "MonkeyC": ("#762a83", "#984ea3", "#c994c7"),
}


def plot_mean_vs_perchannel_comparison(
    all_results: dict[str, dict[str, dict[str, dict]]],
    out_dir: Path,
) -> None:
    """One 2×2 figure: monkeys × areas, mean ARMA vs per-channel summaries."""
    metric_keys = (
        ("rho_avg_test", "mean LFP ← mean MUA"),
        ("rho_ch_mean", "mean(per-ch r)"),
        ("rho_demean_mean", "mean(per-ch r, LFP − avg LFP)"),
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    width = 0.24
    x = np.arange(len(STATE_ORDER))
    subjects = list(STATE_FILES)
    areas = ("PF", "PM")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharey=True)
    for i, subject in enumerate(subjects):
        colors = MONKEY_METRIC_COLORS[subject]
        for j, area in enumerate(areas):
            ax = axes[i, j]
            area_name = "Prefrontal" if area == "PF" else "Premotor"
            for mi, ((key, label), color) in enumerate(zip(metric_keys, colors)):
                offset = (mi - 1) * width
                vals = [
                    all_results[subject][area][state][key] for state in STATE_ORDER
                ]
                ax.bar(
                    x + offset,
                    vals,
                    width * 0.9,
                    color=color,
                    edgecolor="black",
                    linewidth=0.4,
                    label=label,
                )
            ax.set_xticks(x)
            ax.set_xticklabels(STATE_ORDER)
            if j == 0:
                ax.set_ylabel("Test Pearson r")
            ax.set_title(f"{MONKEY_LABELS[subject]} — {area_name}")
            ax.axhline(0, color="gray", linewidth=0.6)
            ax.grid(axis="y", alpha=0.3)
            if i == 0 and j == 0:
                ax.legend(loc="best", fontsize=8)

    fig.tight_layout()
    output = out_dir / "arma_mean_vs_perchannel.png"
    fig.savefig(output, dpi=150, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


def plot_anesthesia_pc1_mean_kernel(
    subject: str,
    area: str,
    kernel: np.ndarray,
    mean_alpha: np.ndarray,
    cfg,
    mea_map: np.ndarray,
    ch_list: list[int],
    all_ch: np.ndarray,
    output: Path,
) -> None:
    """PC1 temporal loading vs mean-mean kernel; spatial PC1 score map."""
    coeff_pc1, ev1, _ = kernel_pca_pc1(kernel)
    out_scores = pc1_scores_per_output(kernel, coeff_pc1)
    lags_s = np.arange(1, KERNEL_LENGTH + 1) * LAG_DT_S

    area_name = "Prefrontal" if area == "PF" else "Premotor"
    color = MONKEY_METRIC_COLORS[subject][1]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    ax.plot(lags_s, coeff_pc1, color="#333333", lw=2, label="PC1 coeff (per-ch kernels)")
    ax.plot(
        lags_s,
        mean_alpha,
        color=color,
        lw=2,
        ls="--",
        label="mean LFP ← mean MUA kernel",
    )
    ax.axhline(0, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("Lag (s)")
    ax.set_ylabel("PC1 loading / kernel weight")
    ax.set_title(f"Temporal PC1 ({ev1:.0%} var.)")
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.3)

    ax = axes[1]
    measure = np.zeros_like(mea_map, dtype=float)
    for ch, k_idx in zip(cfg.ch_out, cfg.kernel_targets):
        r, c = np.where(mea_map == ch)
        if r.size:
            measure[r[0], c[0]] = out_scores[k_idx]
    im = ax.imshow(measure, cmap="RdYlBu_r", origin="upper")
    for ch in ch_list:
        r, c = np.where(mea_map == ch)
        if r.size:
            ax.text(c[0], r[0], str(ch), ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="|PC1 score| per output ch")
    ax.set_title("PC1 score per output channel")
    ax.set_xlabel("X MEA")
    ax.set_ylabel("Y MEA")

    fig.suptitle(
        f"{MONKEY_LABELS[subject]} — {area_name} — {ANESTHESIA_STATE}",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


def run_pc1_plots(scripts_root: Path, out_dir: Path) -> None:
    """Anesthesia PC1 plots for mean-mean kernels (per monkey / area)."""
    global SCRIPTS_ROOT
    SCRIPTS_ROOT = scripts_root
    out_dir.mkdir(parents=True, exist_ok=True)

    for subject in STATE_FILES:
        data = load_subject_data(subject, scripts_root=scripts_root)
        transient = STATE_FILES[subject][ANESTHESIA_STATE]
        for area in ("PF", "PM"):
            cfg = area_config(subject, area)
            mua, lfp = select_mua_lfp(data, cfg)
            spatial_cfg = AREA_CONFIG[subject][area]
            kernels, _ = load_kernel_mat(kernel_path(subject, area, transient))
            mean_alpha = fit_mean_kernel_alpha(mua, lfp, cfg, transient, data["time"])
            all_ch = pf_all_ch() if area == "PF" else pm_all_ch()
            plot_anesthesia_pc1_mean_kernel(
                subject,
                area,
                kernels,
                mean_alpha,
                cfg,
                spatial_cfg["map"],
                spatial_cfg["ch"],
                all_ch,
                out_dir / f"arma_pc1_anesthesia_{subject}_{area}.png",
            )


def results_from_csv(csv_path: Path) -> dict[str, dict[str, dict[str, dict]]]:
    """Rebuild minimal all_results dict from summary CSV (for plot-only runs)."""
    all_results: dict[str, dict[str, dict[str, dict]]] = {
        s: {"PF": {}, "PM": {}} for s in STATE_FILES
    }
    with csv_path.open(newline="") as f:
        for row in csv.DictReader(f):
            all_results[row["subject"]][row["area"]][row["state"]] = {
                k: float(row[k])
                for k in ("rho_avg_test", "rho_ch_mean", "rho_demean_mean")
            }
    return all_results


def plot_boxplot_comparison(
    all_results: dict[str, dict[str, dict[str, dict]]],
    area: str,
    output: Path,
) -> None:
    """Boxplot of per-channel test rho per state, with average markers."""
    fig, ax = plt.subplots(figsize=(9, 5))
    colors_subj = {"MonkeyB": "#4daf4a", "MonkeyC": "#984ea3"}
    width = 0.35

    for si, subject in enumerate(("MonkeyB", "MonkeyC")):
        offset = (si - 0.5) * width
        positions = np.arange(1, len(STATE_ORDER) + 1) + offset
        values = [
            all_results[subject][area][state]["perf_ch_values"] for state in STATE_ORDER
        ]
        bp = ax.boxplot(
            values,
            positions=positions,
            widths=width * 0.85,
            patch_artist=True,
            notch=True,
            showfliers=True,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(colors_subj[subject])
            patch.set_alpha(0.55)
        for pos, state in zip(positions, STATE_ORDER):
            ax.scatter(
                [pos],
                [all_results[subject][area][state]["rho_avg_test"]],
                marker="D",
                s=70,
                color="black",
                zorder=5,
            )

    ax.set_xticks(np.arange(1, len(STATE_ORDER) + 1))
    ax.set_xticklabels(STATE_ORDER)
    ax.set_ylabel("Test Pearson r")
    area_name = "Prefrontal" if area == "PF" else "Premotor"
    ax.set_title(
        f"{area_name} — per-channel ARMA (box) vs mean LFP ← mean MUA (◆)"
    )
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=colors_subj[s], alpha=0.55)
        for s in ("MonkeyB", "MonkeyC")
    ]
    handles.append(
        plt.Line2D(
            [0], [0], marker="D", color="black", linestyle="", markersize=8, label=""
        )
    )
    ax.legend(handles, ["Monkey B", "Monkey C", "avg LFP ← avg MUA"], loc="best")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


def run(scripts_root: Path, out_dir: Path) -> None:
    global SCRIPTS_ROOT
    SCRIPTS_ROOT = scripts_root
    out_dir.mkdir(parents=True, exist_ok=True)

    all_results: dict[str, dict[str, dict[str, dict]]] = {
        s: {"PF": {}, "PM": {}} for s in STATE_FILES
    }
    csv_rows: list[dict] = []

    for subject in STATE_FILES:
        print(f"\n=== {MONKEY_LABELS[subject]} ===", flush=True)
        print("  Loading data for average ARMA...", flush=True)
        data = load_subject_data(subject, scripts_root=scripts_root)

        for area in ("PF", "PM"):
            cfg = area_config(subject, area)
            mua, lfp = select_mua_lfp(data, cfg)
            spatial_cfg = AREA_CONFIG[subject][area]
            mea_map = spatial_cfg["map"]
            ch_list = spatial_cfg["ch"]
            all_ch = pf_all_ch() if area == "PF" else pm_all_ch()

            area_results: dict[str, dict] = {}
            for state, transient in STATE_FILES[subject].items():
                print(f"  {area} / {state} (transient={transient}s)...", flush=True)
                mat_path = kernel_path(subject, area, transient)
                perf_test = load_performance_test(mat_path)
                kernels, _ = load_kernel_mat(mat_path)
                ndx_lp, ndx_tp = test_indices(data["time"], transient)
                in_train = build_lagged_input_fast(mua, cfg.mua_rows, ndx_lp)
                in_test = build_lagged_input_fast(mua, cfg.mua_rows, ndx_tp)
                _, _, rho_avg_te, p_avg_te = average_arma(
                    mua, lfp, cfg, transient, data["time"]
                )
                rho_mean_pred, p_mean_pred = mean_from_per_channel_predictions(
                    kernels, lfp, cfg, in_test, ndx_tp
                )
                perf_demean = per_channel_arma(
                    lfp,
                    cfg,
                    in_train,
                    in_test,
                    ndx_lp,
                    ndx_tp,
                    demean_lfp=True,
                )
                ch_vals = channel_values_from_perf(perf_test, cfg)
                ch_vals_demean = channel_values_from_perf(perf_demean, cfg)
                area_results[state] = {
                    "perf_test": perf_test,
                    "perf_demean": perf_demean,
                    "perf_ch_values": ch_vals,
                    "perf_demean_values": ch_vals_demean,
                    "rho_avg_test": rho_avg_te,
                    "p_avg_test": p_avg_te,
                    "rho_mean_pred": rho_mean_pred,
                    "p_mean_pred": p_mean_pred,
                    "rho_ch_mean": float(np.mean(ch_vals)),
                    "rho_demean_mean": float(np.mean(ch_vals_demean)),
                }
                csv_rows.append(
                    {
                        "subject": subject,
                        "area": area,
                        "state": state,
                        "transient": transient,
                        "rho_avg_test": rho_avg_te,
                        "p_avg_test": p_avg_te,
                        "rho_mean_pred": rho_mean_pred,
                        "p_mean_pred": p_mean_pred,
                        "rho_ch_mean": float(np.mean(ch_vals)),
                        "rho_ch_median": float(np.median(ch_vals)),
                        "rho_ch_std": float(np.std(ch_vals)),
                        "rho_demean_mean": float(np.mean(ch_vals_demean)),
                        "rho_demean_median": float(np.median(ch_vals_demean)),
                        "rho_demean_std": float(np.std(ch_vals_demean)),
                        "n_channels": len(ch_vals),
                    }
                )
                print(
                    f"    per-ch: mean r={ch_vals.mean():.3f}, "
                    f"median={np.median(ch_vals):.3f}  |  "
                    f"mean ARMA: r={rho_avg_te:.3f}  |  "
                    f"mean(pred): r={rho_mean_pred:.3f}  |  "
                    f"demeaned per-ch mean: r={ch_vals_demean.mean():.3f}"
                )

            all_results[subject][area] = area_results
            fig_path = out_dir / f"arma_spatial_{subject}_{area}.png"
            plot_state_row(
                area_results, mea_map, ch_list, all_ch, subject, area, fig_path
            )
            plt.close("all")

    plot_mean_vs_perchannel_comparison(all_results, out_dir)
    run_pc1_plots(scripts_root, out_dir)
    for area in ("PF", "PM"):
        plot_boxplot_comparison(
            all_results,
            area,
            out_dir / f"arma_boxplot_{area}.png",
        )

    csv_path = out_dir / "arma_average_summary.csv"
    fields = list(csv_rows[0].keys())
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(csv_rows)
    print(f"\nSaved summary CSV: {csv_path}")


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
        default=Path(__file__).resolve().parent / "figures",
    )
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Regenerate comparison bar plots from existing arma_average_summary.csv",
    )
    parser.add_argument(
        "--pc1-only",
        action="store_true",
        help="Generate anesthesia PC1 / mean-mean kernel plots only",
    )
    args = parser.parse_args()
    if args.pc1_only:
        run_pc1_plots(args.scripts_root, args.out_dir)
    elif args.plots_only:
        csv_path = args.out_dir / "arma_average_summary.csv"
        plot_mean_vs_perchannel_comparison(results_from_csv(csv_path), args.out_dir)
    else:
        run(args.scripts_root, args.out_dir)
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
