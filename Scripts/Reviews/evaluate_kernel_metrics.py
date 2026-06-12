#!/usr/bin/env python3
"""Evaluate precomputed ARMA kernels on held-out test windows.

Loads kernels from Scripts/MonkeyB|MonkeyC/Data/Kernel, applies them to the test period defined
in each ARMA_*.m script, and reports correlation, RMSE, R², and comparison
against a mean (null) predictor.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import scipy.io
from scipy import stats

SCRIPTS_ROOT = Path(__file__).resolve().parent.parent

KERNEL_LENGTH = 69
LEARNING_PERIOD = 200
TEST_PERIOD = 50
SMOOTHING_WINDOW = 1

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

AREAS = ("PF", "PM", "PF2PM", "PM2PF")


@dataclass(frozen=True)
class AreaConfig:
    """Channel layout and MUA/LFP indexing for one ARMA script."""

    ch_in: np.ndarray  # 1-based channel IDs used as input
    ch_out: np.ndarray  # 1-based channel IDs used as output / kernel target index
    mua_source: str  # 'pf' | 'pm' | 'by_channel' — which MUA buffer ARMA_*.m uses
    mua_rows: np.ndarray  # 0-based rows into that buffer (one per ch_in)
    lfp_source: str  # 'pf' | 'pm'
    lfp_rows: np.ndarray  # 0-based rows into LFP buffer (one per ch_out)
    kernel_targets: np.ndarray  # 0-based index into Kernel[:, :, target]


def area_config(subject: str, area: str) -> AreaConfig:
    """Mirror channel choices in Scripts/{subject}/ARMA_*.m."""
    if subject == "MonkeyC":
        if area == "PF":
            ch = np.setdiff1d(np.arange(1, 97), [50, 62])
            return AreaConfig(
                ch_in=ch,
                ch_out=ch,
                mua_source="pf",
                mua_rows=ch - 1,
                lfp_source="pf",
                lfp_rows=ch - 1,
                kernel_targets=ch - 1,
            )
        if area == "PM":
            ch = np.arange(97, 193)
            return AreaConfig(
                ch_in=ch,
                ch_out=ch,
                mua_source="pm",
                mua_rows=ch - 97,
                lfp_source="pm",
                lfp_rows=ch - 97,
                kernel_targets=ch - 97,
            )
        if area == "PF2PM":
            ch_in = np.setdiff1d(np.arange(1, 97), [50, 62])
            ch_out = np.arange(97, 193)
            return AreaConfig(
                ch_in=ch_in,
                ch_out=ch_out,
                mua_source="by_channel",
                mua_rows=ch_in - 1,
                lfp_source="pm",
                lfp_rows=ch_out - 97,
                kernel_targets=ch_out - 97,
            )
        if area == "PM2PF":
            ch_in = np.arange(97, 193)
            ch_out = np.setdiff1d(np.arange(1, 97), [62, 50])
            return AreaConfig(
                ch_in=ch_in,
                ch_out=ch_out,
                mua_source="by_channel",
                mua_rows=ch_in - 1,
                lfp_source="pf",
                lfp_rows=ch_out - 1,
                kernel_targets=ch_out - 1,
            )
    if subject == "MonkeyB":
        if area == "PF":
            ch = np.setdiff1d(np.arange(1, 97), [66, 20, 88])
            return AreaConfig(
                ch_in=ch,
                ch_out=ch,
                mua_source="pf",
                mua_rows=ch - 1,
                lfp_source="pf",
                lfp_rows=ch - 1,
                kernel_targets=ch - 1,
            )
        if area == "PM":
            ch = np.setdiff1d(np.arange(97, 193), [174, 113, 131, 133])
            return AreaConfig(
                ch_in=ch,
                ch_out=ch,
                mua_source="pm",
                mua_rows=ch - 97,
                lfp_source="pm",
                lfp_rows=ch - 97,
                kernel_targets=ch - 97,
            )
        if area == "PF2PM":
            ch_in = np.setdiff1d(np.arange(1, 97), [66, 20, 88])
            ch_out = np.setdiff1d(np.arange(97, 193), [174, 113, 131, 133])
            return AreaConfig(
                ch_in=ch_in,
                ch_out=ch_out,
                mua_source="by_channel",
                mua_rows=ch_in - 1,
                lfp_source="pm",
                lfp_rows=ch_out - 97,
                kernel_targets=ch_out - 97,
            )
        if area == "PM2PF":
            ch_in = np.setdiff1d(np.arange(97, 193), [174, 113, 131, 133])
            ch_out = np.setdiff1d(np.arange(1, 97), [66, 20, 88])
            return AreaConfig(
                ch_in=ch_in,
                ch_out=ch_out,
                mua_source="by_channel",
                mua_rows=ch_in - 1,
                lfp_source="pf",
                lfp_rows=ch_out - 1,
                kernel_targets=ch_out - 1,
            )
    raise ValueError(f"Unknown subject/area: {subject}/{area}")


def load_subject_data(subject: str, scripts_root: Path = SCRIPTS_ROOT) -> dict[str, np.ndarray]:
    """Load log-MUA and LFP buffers matching ARMA_*.m layouts."""
    data_dir = scripts_root / subject
    pre = scipy.io.loadmat(data_dir / "PreProcData.mat")
    log_shift = pre["DataSet"]["LogMUAshift"][0, 0].ravel()

    mat_path = data_dir / "MEAMUALFP.mat"
    with h5py.File(mat_path, "r") as f:
        time = np.asarray(f["MEAMUA"]["time"]).squeeze()
        mua_raw = np.asarray(f["MEAMUA"]["values"]).T
        lfp_raw = np.asarray(f["MEALFP"]["values"]).T

    n_ch = mua_raw.shape[0]
    mua_log = np.empty_like(mua_raw, dtype=float)
    for c in range(n_ch):
        mua_log[c] = np.log(mua_raw[c]) - log_shift[c]

    # PF scripts: MUA(c,:) stores channel c for c=1..96
    mua_pf = mua_log[:96]
    # PM scripts: MUA(c,:) stores channel 96+c (compact 97..192)
    mua_pm = mua_log[96:192]
    # Inter-area scripts: MUA(channel_id,:) indexed by global channel number
    mua_by_channel = mua_log

    lfp_pf = lfp_raw[:96]
    lfp_pm = lfp_raw[96:192]

    return {
        "time": time,
        "mua_pf": mua_pf,
        "mua_pm": mua_pm,
        "mua_by_channel": mua_by_channel,
        "lfp_pf": lfp_pf,
        "lfp_pm": lfp_pm,
    }


def select_mua_lfp(data: dict[str, np.ndarray], cfg: AreaConfig) -> tuple[np.ndarray, np.ndarray]:
    mua_map = {
        "pf": data["mua_pf"],
        "pm": data["mua_pm"],
        "by_channel": data["mua_by_channel"],
    }
    lfp_map = {"pf": data["lfp_pf"], "pm": data["lfp_pm"]}
    return mua_map[cfg.mua_source], lfp_map[cfg.lfp_source]


def test_indices(time: np.ndarray, transient_period: float) -> tuple[np.ndarray, np.ndarray]:
    """Learning and test sample indices (MATLAB ARMA scripts)."""
    t0 = time[0]
    ndx_lp = np.where(
        (time >= t0 + transient_period)
        & (time <= t0 + transient_period + LEARNING_PERIOD)
    )[0]
    ndx_tp = np.where(
        (time >= t0 + transient_period + LEARNING_PERIOD)
        & (time <= t0 + transient_period + LEARNING_PERIOD + TEST_PERIOD)
    )[0]
    return ndx_lp, ndx_tp


def build_lagged_input(
    mua: np.ndarray, row_idx: np.ndarray, sample_idx: np.ndarray
) -> np.ndarray:
    """Stack lagged MUA blocks: shape (n_ch * kernel_length, n_samples)."""
    offsets = np.arange(1, KERNEL_LENGTH + 1, SMOOTHING_WINDOW)
    blocks = [mua[row_idx][:, sample_idx - off] for off in offsets]
    return np.vstack(blocks)


def kernel_to_alpha(kernel_slice: np.ndarray) -> np.ndarray:
    """Flatten KernelLength x n_ch matrix to Alpha row vector (MATLAB order)."""
    return kernel_slice.T.reshape(-1, order="F")


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    """Pearson r, RMSE, R², and null (mean) baseline on the same window."""
    y_true = np.asarray(y_true, dtype=float).ravel()
    y_pred = np.asarray(y_pred, dtype=float).ravel()
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true = y_true[mask]
    y_pred = y_pred[mask]
    if y_true.size < 3:
        return {
            "rho": np.nan,
            "p_rho": np.nan,
            "rmse": np.nan,
            "r2": np.nan,
            "rmse_null": np.nan,
            "r2_null": 0.0,
            "delta_r2": np.nan,
            "rmse_gain": np.nan,
        }

    rho, p_rho = stats.pearsonr(y_true, y_pred)
    resid = y_true - y_pred
    rmse = float(np.sqrt(np.mean(resid**2)))

    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    ss_res = float(np.sum(resid**2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan

    y_null = np.full_like(y_true, y_true.mean())
    rmse_null = float(np.sqrt(np.mean((y_true - y_null) ** 2)))
    r2_null = 0.0
    delta_r2 = r2 - r2_null
    rmse_gain = (rmse_null - rmse) / rmse_null if rmse_null > 0 else np.nan

    return {
        "rho": float(rho),
        "p_rho": float(p_rho),
        "rmse": rmse,
        "r2": r2,
        "rmse_null": rmse_null,
        "r2_null": r2_null,
        "delta_r2": float(delta_r2),
        "rmse_gain": float(rmse_gain),
    }


def kernel_mat_path(
    subject: str, area: str, transient: int, scripts_root: Path = SCRIPTS_ROOT
) -> Path:
    return (
        scripts_root
        / subject
        / "Data"
        / "Kernel"
        / f"KernelAndPerf{area}_trans{transient}.mat"
    )


def load_kernel_mat(mat_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load only Kernel and PerformanceTest (skip large reconstruction arrays)."""
    mat = scipy.io.loadmat(
        mat_path, variable_names=["Kernel", "PerformanceTest"], simplify_cells=True
    )
    kernels = mat["Kernel"]
    perf = np.asarray(mat["PerformanceTest"]).ravel()
    return kernels, perf


def evaluate_condition(
    subject: str,
    area: str,
    state: str,
    transient: int,
    data: dict[str, np.ndarray],
    cfg: AreaConfig,
    scripts_root: Path = SCRIPTS_ROOT,
) -> tuple[list[dict], dict]:
    """Apply stored kernels on the test window for all output channels."""
    mat_path = kernel_mat_path(subject, area, transient, scripts_root)
    kernels, perf_stored = load_kernel_mat(mat_path)

    mua, lfp = select_mua_lfp(data, cfg)
    _, ndx_tp = test_indices(data["time"], transient)
    input_test = build_lagged_input(mua, cfg.mua_rows, ndx_tp)

    rows: list[dict] = []
    for out_ch, k_idx, lfp_row in zip(cfg.ch_out, cfg.kernel_targets, cfg.lfp_rows):
        y_true = lfp[lfp_row, ndx_tp + 1]

        alpha = kernel_to_alpha(kernels[:, :, k_idx])
        y_pred = alpha @ input_test

        m = regression_metrics(y_true, y_pred)
        rows.append(
            {
                "subject": subject,
                "area": area,
                "state": state,
                "transient": transient,
                "channel": int(out_ch),
                "rho": m["rho"],
                "p_rho": m["p_rho"],
                "rho_stored": float(perf_stored[k_idx]),
                "rho_diff": m["rho"] - float(perf_stored[k_idx]),
                "rmse": m["rmse"],
                "r2": m["r2"],
                "rmse_null": m["rmse_null"],
                "r2_null": m["r2_null"],
                "delta_r2": m["delta_r2"],
                "rmse_gain": m["rmse_gain"],
                "beats_null_r2": m["delta_r2"] > 0,
                "sig_rho_05": m["p_rho"] < 0.05,
            }
        )

    rhos = np.array([r["rho"] for r in rows])
    r2s = np.array([r["r2"] for r in rows])
    summary = {
        "subject": subject,
        "area": area,
        "state": state,
        "transient": transient,
        "n_channels": len(rows),
        "rho_mean": float(np.nanmean(rhos)),
        "rho_median": float(np.nanmedian(rhos)),
        "r2_mean": float(np.nanmean(r2s)),
        "r2_median": float(np.nanmedian(r2s)),
        "rmse_mean": float(np.nanmean([r["rmse"] for r in rows])),
        "frac_sig_rho": float(np.mean([r["sig_rho_05"] for r in rows])),
        "frac_beats_null": float(np.mean([r["beats_null_r2"] for r in rows])),
        "max_abs_rho_diff": float(np.nanmax(np.abs([r["rho_diff"] for r in rows]))),
    }

    if len(rhos) > 0:
        wstat, p_wilcox = stats.wilcoxon(rhos, alternative="greater")
        tstat, p_t = stats.ttest_1samp(rhos, 0.0, alternative="greater")
        summary["wilcoxon_rho_vs_0_stat"] = float(wstat)
        summary["wilcoxon_rho_vs_0_p"] = float(p_wilcox)
        summary["ttest_rho_vs_0_stat"] = float(tstat)
        summary["ttest_rho_vs_0_p"] = float(p_t)
        wstat_r2, p_wilcox_r2 = stats.wilcoxon(r2s, alternative="greater")
        summary["wilcoxon_r2_vs_null_stat"] = float(wstat_r2)
        summary["wilcoxon_r2_vs_null_p"] = float(p_wilcox_r2)

    return rows, summary


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run(scripts_root: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    all_channel_rows: list[dict] = []
    all_summary_rows: list[dict] = []

    for subject in STATE_FILES:
        print(f"Loading {subject} data...")
        data = load_subject_data(subject, scripts_root)

        for area in AREAS:
            cfg = area_config(subject, area)
            for state, transient in STATE_FILES[subject].items():
                print(f"  {area} / {state} (trans{transient})...")
                ch_rows, summary = evaluate_condition(
                    subject, area, state, transient, data, cfg, scripts_root
                )
                all_channel_rows.extend(ch_rows)
                all_summary_rows.append(summary)

    channel_fields = [
        "subject",
        "area",
        "state",
        "transient",
        "channel",
        "rho",
        "p_rho",
        "rho_stored",
        "rho_diff",
        "rmse",
        "r2",
        "rmse_null",
        "r2_null",
        "delta_r2",
        "rmse_gain",
        "beats_null_r2",
        "sig_rho_05",
    ]
    summary_fields = [
        "subject",
        "area",
        "state",
        "transient",
        "n_channels",
        "rho_mean",
        "rho_median",
        "r2_mean",
        "r2_median",
        "rmse_mean",
        "frac_sig_rho",
        "frac_beats_null",
        "max_abs_rho_diff",
        "wilcoxon_rho_vs_0_stat",
        "wilcoxon_rho_vs_0_p",
        "ttest_rho_vs_0_stat",
        "ttest_rho_vs_0_p",
        "wilcoxon_r2_vs_null_stat",
        "wilcoxon_r2_vs_null_p",
    ]

    ch_path = out_dir / "kernel_metrics_per_channel.csv"
    sum_path = out_dir / "kernel_metrics_summary.csv"
    write_csv(ch_path, all_channel_rows, channel_fields)
    write_csv(sum_path, all_summary_rows, summary_fields)

    print(f"\nSaved per-channel metrics: {ch_path}")
    print(f"Saved summary metrics:     {sum_path}")
    print_summary_table(all_summary_rows)


def print_summary_table(rows: list[dict]) -> None:
    print("\n" + "=" * 100)
    print(
        "Summary (test window): mean rho | mean R² | frac sig rho | frac beats null | "
        "Wilcoxon p(rho>0)"
    )
    print("=" * 100)
    for r in rows:
        p_rho = r.get("wilcoxon_rho_vs_0_p", np.nan)
        p_str = f"{p_rho:.2e}" if p_rho < 1e-3 else f"{p_rho:.4g}"
        print(
            f"{r['subject']:9s} {r['area']:6s} {r['state']:12s} "
            f"rho={r['rho_mean']:.3f}  R2={r['r2_mean']:.3f}  "
            f"sig={r['frac_sig_rho']:.2f}  beat_null={r['frac_beats_null']:.2f}  "
            f"p_rho={p_str}  max|rho_diff|={r['max_abs_rho_diff']:.4g}"
        )


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
        default=Path(__file__).resolve().parent / "kernel_metrics",
    )
    args = parser.parse_args()
    run(args.scripts_root, args.out_dir)


if __name__ == "__main__":
    main()
