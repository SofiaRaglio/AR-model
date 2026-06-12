#!/usr/bin/env python3
"""Python port of detectingUpAndDownStates.m — up/down state durations from MEAMUALFP."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from pathlib import Path

import h5py
import numpy as np
from scipy import stats
from scipy.signal import butter, filtfilt
from sklearn.decomposition import PCA

SCRIPTS_ROOT = Path(__file__).resolve().parent.parent

DEFAULT_DATASETS = {
    "MonkeyB": SCRIPTS_ROOT / "MonkeyB" / "MEAMUALFP.mat",
    "MonkeyC": SCRIPTS_ROOT / "MonkeyC" / "MEAMUALFP.mat",
}

SUBJECT_LABELS = {
    "MonkeyB": "Monkey B",
    "MonkeyC": "Monkey C",
}

# [0 250] + transient: Monkey B uses 3000 s, Monkey C uses 4000 s
DEFAULT_TIME_RANGES: dict[str, tuple[float, float]] = {
    "MonkeyB": (3000.0, 3250.0),
    "MonkeyC": (4000.0, 4250.0),
}

# Defaults match detectingUpAndDownStates.m (TimeRange, Threshold, lpCutOffFreq, etc.)
DEFAULT_PARAMS = None  # set after DetectionParams is defined


@dataclass
class DetectionParams:
    """Tunable detection / duration parameters."""

    time_start: float = 3000.0
    time_end: float = 3250.0
    up_threshold: float = 0.4  # MATLAB: UpDownDect > 0.4
    min_duration: float = 0.3  # MATLAB: Threshold for kept epochs
    lp_cutoff_hz: float = 12.5
    filter_order: int = 4

    @property
    def time_range(self) -> tuple[float, float]:
        return (self.time_start, self.time_end)


DEFAULT_PARAMS = DetectionParams()


def add_detection_args(parser: argparse.ArgumentParser) -> None:
    g = parser.add_argument_group("detection parameters (MATLAB defaults in brackets)")
    g.add_argument(
        "--time-start",
        type=float,
        default=None,
        help="Window start (s); default: Monkey B 3000, Monkey C 4000",
    )
    g.add_argument(
        "--time-end",
        type=float,
        default=None,
        help="Window end (s); default: start + 250 s",
    )
    g.add_argument(
        "--up-threshold",
        type=float,
        default=0.4,
        help="Up/down boundary on normalized detector [0.4]",
    )
    g.add_argument(
        "--min-duration",
        type=float,
        default=0.3,
        help="Drop epochs shorter than this (s) [0.3]",
    )
    g.add_argument("--lp-cutoff", type=float, default=12.5, help="MUA low-pass (Hz) [12.5]")
    g.add_argument("--filter-order", type=int, default=4, help="Butterworth order [4]")


def detection_params_from_args(args: argparse.Namespace) -> DetectionParams:
    return DetectionParams(
        time_start=args.time_start if args.time_start is not None else 3000.0,
        time_end=args.time_end if args.time_end is not None else 3250.0,
        up_threshold=args.up_threshold,
        min_duration=args.min_duration,
        lp_cutoff_hz=args.lp_cutoff,
        filter_order=args.filter_order,
    )


def time_range_overridden_from_cli(args: argparse.Namespace) -> bool:
    return args.time_start is not None or args.time_end is not None


def params_for_mat(
    mat_path: Path,
    params: DetectionParams,
    cli_time_override: bool,
) -> DetectionParams:
    """Apply per-subject default time window unless CLI overrides it."""
    if cli_time_override:
        return params
    subject = mat_path.parent.name
    if subject in DEFAULT_TIME_RANGES:
        t0, t1 = DEFAULT_TIME_RANGES[subject]
        return replace(params, time_start=t0, time_end=t1)
    return params


def load_meamualfp_struct(
    mat_path: Path, name: str, row_idx: np.ndarray | None = None
) -> dict[str, np.ndarray]:
    """Load MEAMUA or MEALFP from a MATLAB v7.3 (HDF5) file."""
    with h5py.File(mat_path, "r") as f:
        g = f[name]
        time_ds = g["time"]
        if row_idx is None:
            time = np.asarray(time_ds).squeeze()
            values = np.asarray(g["values"])
        else:
            time = np.asarray(time_ds[row_idx]).squeeze()
            values = np.asarray(g["values"][row_idx, :])
        dt = float(np.asarray(g["dt"]).squeeze())
    return {"time": time, "values": values.T, "dt": dt}  # channels x time


def lowpass_rows(x: np.ndarray, cutoff_hz: float, fs_hz: float, order: int = 4) -> np.ndarray:
    """Zero-phase lowpass along the time axis (MATLAB lowpass-style)."""
    wn = cutoff_hz / (fs_hz / 2.0)
    b, a = butter(order, wn, btype="low")
    return filtfilt(b, a, x, axis=1)


def pca_first_component(x_time_channels: np.ndarray) -> np.ndarray:
    """First PCA score; x has shape (n_time, n_features)."""
    pca = PCA(n_components=1)
    return pca.fit_transform(x_time_channels).ravel()


def align_detector(detector: np.ndarray, reference: np.ndarray) -> np.ndarray:
    detector = detector / np.std(detector)
    reference = reference / np.std(reference)
    r = np.corrcoef(detector, reference)[0, 1]
    return detector * np.sign(r)


def up_state_detector(
    data_channels_time: np.ndarray, reference: np.ndarray | None = None
) -> np.ndarray:
    """PCA-based up/down detector; data shape (n_channels, n_time)."""
    detector = pca_first_component(data_channels_time.T)
    ref = data_channels_time.mean(axis=0) if reference is None else reference
    return align_detector(detector, ref)


def state_durations(
    up_down: np.ndarray,
    dt: float,
    threshold: float = 0.4,
    min_duration: float = 0.3,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    transitions = (up_down > threshold).astype(np.int8)
    ndx = np.where(np.diff(transitions) != 0)[0] + 1
    if ndx.size == 0:
        return np.array([]), np.array([]), {
            "mean_up": np.nan,
            "std_up": np.nan,
            "mean_down": np.nan,
            "std_down": np.nan,
            "n_up": 0,
            "n_down": 0,
        }

    # duration[i] is the interval starting at transition ndx[i] (MATLAB indexing)
    duration = np.diff(ndx) * dt
    val = transitions[ndx[:-1]]
    keep = duration >= min_duration
    duration = duration[keep]
    val = val[keep]

    up_dur = duration[val == 1]
    down_dur = duration[val == 0]

    stats = {
        "mean_up": float(np.mean(up_dur)) if up_dur.size else np.nan,
        "std_up": float(np.std(up_dur)) if up_dur.size else np.nan,
        "mean_down": float(np.mean(down_dur)) if down_dur.size else np.nan,
        "std_down": float(np.std(down_dur)) if down_dur.size else np.nan,
        "n_up": int(up_dur.size),
        "n_down": int(down_dur.size),
    }
    return up_dur, down_dur, stats


def build_detector(
    modality: str,
    lp_z: np.ndarray,
    lfp_z: np.ndarray,
) -> np.ndarray:
    """
    Build up/down trace.

    - mua: PCA on low-pass log-MUA only
    - lfp: PCA on LFP only
    - combined: PCA on [lp MUA; LFP] — same as detectingUpAndDownStates.m
    """
    if modality == "mua":
        return up_state_detector(lp_z)
    if modality == "lfp":
        return up_state_detector(lfp_z)
    if modality == "combined":
        detector = pca_first_component(np.vstack([lp_z, lfp_z]).T)
        return align_detector(detector, lp_z.mean(axis=0))
    raise ValueError(f"Unknown modality: {modality}")


def run_modality_up_durations(
    mat_path: Path,
    params: DetectionParams | None = None,
    modalities: tuple[str, ...] = ("mua", "lfp", "combined"),
    cli_time_override: bool = False,
) -> dict:
    """Up-state durations per detector type (PF and PM)."""
    params = params_for_mat(mat_path, params or DEFAULT_PARAMS, cli_time_override)
    with h5py.File(mat_path, "r") as f:
        mua_time = np.asarray(f["MEAMUA"]["time"]).squeeze()
        dt = float(np.asarray(f["MEAMUA"]["dt"]).squeeze())

    ndx = np.where(
        (mua_time >= params.time_start) & (mua_time <= params.time_end)
    )[0]
    meamua = load_meamualfp_struct(mat_path, "MEAMUA", ndx)
    mealfp = load_meamualfp_struct(mat_path, "MEALFP", ndx)
    fs = 1.0 / dt

    pf_channels = np.setdiff1d(np.arange(1, 97), [20, 66, 88]) - 1
    pm_channels = np.setdiff1d(np.arange(97, 193), [174, 113, 131, 133]) - 1

    out: dict = {"dt": dt, "subject": mat_path.parent.name, "params": params}
    for region, ch in (("pf", pf_channels), ("pm", pm_channels)):
        z = np.log(meamua["values"][ch])
        z = z - z.mean(axis=1, keepdims=True)
        lp_z = lowpass_rows(
            z, params.lp_cutoff_hz, fs, order=params.filter_order
        )
        lfp_z = mealfp["values"][ch]

        for modality in modalities:
            detector = build_detector(modality, lp_z, lfp_z)
            up_dur, _, up_stats = state_durations(
                detector,
                dt,
                threshold=params.up_threshold,
                min_duration=params.min_duration,
            )
            key = f"{region}_{modality}"
            out[f"{key}_up_durations"] = up_dur
            out[f"{key}_stats"] = up_stats

    return out


def compare_areas_test(
    pf_durations: np.ndarray, pm_durations: np.ndarray
) -> dict[str, float | str]:
    """PF vs PM comparison on up-state duration samples."""
    if pf_durations.size < 2 or pm_durations.size < 2:
        return {"test": "insufficient data", "n_pf": pf_durations.size, "n_pm": pm_durations.size}

    t_stat, t_p = stats.ttest_ind(pf_durations, pm_durations, equal_var=False)
    u_stat, u_p = stats.mannwhitneyu(
        pf_durations, pm_durations, alternative="two-sided"
    )
    return {
        "test": "Welch t-test; Mann-Whitney U",
        "n_pf": int(pf_durations.size),
        "n_pm": int(pm_durations.size),
        "pf_mean": float(np.mean(pf_durations)),
        "pm_mean": float(np.mean(pm_durations)),
        "t_stat": float(t_stat),
        "t_p": float(t_p),
        "u_stat": float(u_stat),
        "u_p": float(u_p),
    }


def plot_modality_comparison(
    results: dict[str, dict],
    out_path: Path,
    modalities: tuple[str, ...] = ("mua", "lfp", "combined"),
) -> None:
    """Bar plot of mean up-state duration: modality × PF/PM × subjects."""
    import matplotlib.pyplot as plt

    modality_labels = {
        "mua": "MUA only",
        "lfp": "LFP only",
        "combined": "MUA+LFP\n(MATLAB)",
    }

    subjects = list(results.keys())
    modalities = tuple(m for m in modalities if f"pf_{m}_stats" in results[subjects[0]])
    regions = ("pf", "pm")

    n_groups = len(subjects) * len(modalities)
    x = np.arange(n_groups)
    width = 0.35
    colors = {"pf": "#2166ac", "pm": "#b2182b"}

    fig, ax = plt.subplots(figsize=(max(10, n_groups * 1.2), 5))
    labels: list[str] = []
    bar_tops: list[float] = []

    for gi, subject in enumerate(subjects):
        for mi, modality in enumerate(modalities):
            pos = gi * len(modalities) + mi
            label = SUBJECT_LABELS.get(subject, subject)
            labels.append(f"{label}\n{modality_labels.get(modality, modality)}")
            group_heights: list[float] = []
            for ri, region in enumerate(regions):
                key = f"{region}_{modality}"
                s = results[subject][f"{key}_stats"]
                mean = s["mean_up"]
                n = s["n_up"]
                sem = s["std_up"] / np.sqrt(n) if n else np.nan
                offset = (ri - 0.5) * width
                ax.bar(
                    pos + offset,
                    mean,
                    width,
                    yerr=sem,
                    capsize=3,
                    color=colors[region],
                    edgecolor="black",
                    linewidth=0.5,
                )
                top = mean + sem if np.isfinite(sem) else mean
                group_heights.append(top)

            pf_dur = results[subject][f"pf_{modality}_up_durations"]
            pm_dur = results[subject][f"pm_{modality}_up_durations"]
            test = compare_areas_test(pf_dur, pm_dur)
            p_val = test.get("u_p", np.nan)
            if np.isfinite(p_val):
                y_line = max(group_heights) * 1.05 + 0.02
                ax.plot(
                    [pos - width / 2, pos + width / 2],
                    [y_line, y_line],
                    "k-",
                    linewidth=0.8,
                )
                sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
                ax.text(pos, y_line + 0.01, sig, ha="center", va="bottom", fontsize=9)
                bar_tops.append(y_line)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Up-state duration (s)")
    ax.set_title("Up-state duration by modality, area, and subject")
    if bar_tops:
        ax.set_ylim(0, max(bar_tops) * 1.15)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=colors["pf"]),
        plt.Rectangle((0, 0), 1, 1, color=colors["pm"]),
    ]
    ax.legend(handles, ["PF", "PM"], loc="upper right")
    ax.text(
        0.01,
        0.99,
        "Brackets: Mann-Whitney PF vs PM (* p<0.05)",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
    )
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def print_detection_params(params: DetectionParams) -> None:
    print("Detection parameters:")
    print(f"  time window: {params.time_start} – {params.time_end} s")
    print(f"  up_threshold: {params.up_threshold}")
    print(f"  min_duration: {params.min_duration} s")
    print(f"  lp_cutoff: {params.lp_cutoff_hz} Hz, filter order: {params.filter_order}")
    print(
        "  Time defaults: Monkey B 3000–3250 s, Monkey C 4000–4250 s "
        "(override with --time-start / --time-end)."
    )


def print_area_comparison_stats(
    results: dict[str, dict],
    modalities: tuple[str, ...] = ("mua", "lfp", "combined"),
) -> None:
    """Print PF vs PM tests for each subject and modality."""
    if results:
        p0 = next(iter(results.values())).get("params")
        if isinstance(p0, DetectionParams):
            print_detection_params(p0)
    for subject in results:
        label = SUBJECT_LABELS.get(subject, subject)
        print(f"\n=== {label}: PF vs PM (up-state duration) ===")
        for modality in modalities:
            if f"pf_{modality}_up_durations" not in results[subject]:
                continue
            pf = results[subject]["pf_" + modality + "_up_durations"]
            pm = results[subject]["pm_" + modality + "_up_durations"]
            res = compare_areas_test(pf, pm)
            print(f"\n  {modality.upper()}:")
            if res.get("test") == "insufficient data":
                print(f"    {res}")
                continue
            print(
                f"    PF mean = {res['pf_mean']:.4f} s (n={res['n_pf']}), "
                f"PM mean = {res['pm_mean']:.4f} s (n={res['n_pm']})"
            )
            print(
                f"    Welch t-test: t = {res['t_stat']:.3f}, p = {res['t_p']:.4g}"
            )
            print(
                f"    Mann-Whitney U: U = {res['u_stat']:.3f}, p = {res['u_p']:.4g}"
            )


def run_area_comparison(
    datasets: dict[str, Path] | None = None,
    out_dir: Path | None = None,
    params: DetectionParams | None = None,
    modalities: tuple[str, ...] = ("mua", "lfp", "combined"),
    cli_time_override: bool = False,
) -> dict[str, dict]:
    datasets = datasets or DEFAULT_DATASETS
    params = params or DEFAULT_PARAMS
    results: dict[str, dict] = {}
    for name, mat_path in datasets.items():
        if not mat_path.is_file():
            raise SystemExit(f"MAT file not found: {mat_path}")
        results[name] = run_modality_up_durations(
            mat_path, params, modalities, cli_time_override
        )

    out_dir = out_dir or Path(__file__).resolve().parent
    plot_path = out_dir / "up_state_duration_mua_lfp_comparison.png"
    plot_modality_comparison(results, plot_path, modalities)
    print_area_comparison_stats(results, modalities)
    print(f"\nSaved figure: {plot_path}")
    return results


def run(
    mat_path: Path,
    params: DetectionParams | None = None,
    cli_time_override: bool = False,
) -> dict:
    params = params_for_mat(mat_path, params or DEFAULT_PARAMS, cli_time_override)
    time_range = params.time_range
    with h5py.File(mat_path, "r") as f:
        mua_time = np.asarray(f["MEAMUA"]["time"]).squeeze()
        dt = float(np.asarray(f["MEAMUA"]["dt"]).squeeze())

    ndx = np.where((mua_time >= time_range[0]) & (mua_time <= time_range[1]))[0]
    meamua = load_meamualfp_struct(mat_path, "MEAMUA", ndx)
    mealfp = load_meamualfp_struct(mat_path, "MEALFP", ndx)
    fs = 1.0 / dt

    pf_channels = np.setdiff1d(np.arange(1, 97), [20, 66, 88]) - 1  # 0-based
    pm_channels = np.setdiff1d(np.arange(97, 193), [174, 113, 131, 133]) - 1

    t = meamua["time"] - meamua["time"][0]

    zpf = np.log(meamua["values"][pf_channels])
    zpm = np.log(meamua["values"][pm_channels])
    zpf = zpf - zpf.mean(axis=1, keepdims=True)
    zpm = zpm - zpm.mean(axis=1, keepdims=True)

    lp_zpf = lowpass_rows(zpf, params.lp_cutoff_hz, fs, order=params.filter_order)
    lp_zpm = lowpass_rows(zpm, params.lp_cutoff_hz, fs, order=params.filter_order)

    lfp_zpf = mealfp["values"][pf_channels]
    lfp_zpm = mealfp["values"][pm_channels]

    lp_lfp_pf = lowpass_rows(lfp_zpf, params.lp_cutoff_hz, fs, order=params.filter_order)
    lp_lfp_pm = lowpass_rows(lfp_zpm, params.lp_cutoff_hz, fs, order=params.filter_order)

    # PF cortex
    pf_up_down = pca_first_component(
        np.vstack([lp_zpf, lfp_zpf]).T
    )
    pf_mean_log_mua = lp_zpf.mean(axis=0)
    pf_up_down = align_detector(pf_up_down, pf_mean_log_mua)

    # PM cortex
    pm_up_down = pca_first_component(
        np.vstack([lp_zpm, lfp_zpm]).T
    )
    pm_mean_log_mua = lp_zpm.mean(axis=0)
    pm_up_down = align_detector(pm_up_down, pm_mean_log_mua)

    pf_up, pf_down, pf_stats = state_durations(
        pf_up_down, dt, params.up_threshold, params.min_duration
    )
    pm_up, pm_down, pm_stats = state_durations(
        pm_up_down, dt, params.up_threshold, params.min_duration
    )

    return {
        "time_s": t,
        "dt": dt,
        "params": params,
        "pf": pf_stats,
        "pm": pm_stats,
        "pf_up_durations": pf_up,
        "pf_down_durations": pf_down,
        "pm_up_durations": pm_up,
        "pm_down_durations": pm_down,
        "lp_zpf": lp_zpf,
        "lp_zpm": lp_zpm,
        "lp_lfp_pf": lp_lfp_pf,
        "lp_lfp_pm": lp_lfp_pm,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mat",
        type=Path,
        default=DEFAULT_DATASETS["MonkeyB"],
        help="Path to MEAMUALFP.mat",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Save summary figures (requires matplotlib)",
    )
    parser.add_argument(
        "--compare-areas",
        action="store_true",
        help=(
            "Plot MUA/LFP up-state durations (PF vs PM) for Monkey B and Monkey C "
            "and run PF–PM statistical tests"
        ),
    )
    parser.add_argument(
        "--compare-out-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for the multi-subject comparison figure",
    )
    parser.add_argument(
        "--compare-modalities",
        nargs="+",
        default=["mua", "lfp", "combined"],
        choices=["mua", "lfp", "combined"],
        help="Detector types in comparison plot (default: all three)",
    )
    add_detection_args(parser)
    args = parser.parse_args()
    params = detection_params_from_args(args)
    cli_time_override = time_range_overridden_from_cli(args)

    if args.compare_areas:
        run_area_comparison(
            out_dir=args.compare_out_dir,
            params=params,
            modalities=tuple(args.compare_modalities),
            cli_time_override=cli_time_override,
        )
        return

    if not args.mat.is_file():
        raise SystemExit(f"MAT file not found: {args.mat}")

    out = run(args.mat, params, cli_time_override)

    print(f"Data: {args.mat}")
    print(f"dt = {out['dt']:.6f} s, window = {out['time_s'][0]:.1f}–{out['time_s'][-1]:.1f} s rel.")
    print()
    for region, label in (("pf", "PF"), ("pm", "PM")):
        s = out[region]
        print(f"{label} cortex:")
        print(
            f"  Up-state:   mean = {s['mean_up']:.4f} s, "
            f"std = {s['std_up']:.4f} s (n = {s['n_up']})"
        )
        print(
            f"  Down-state: mean = {s['mean_down']:.4f} s, "
            f"std = {s['std_down']:.4f} s (n = {s['n_down']})"
        )
        se_up = s["std_up"] / np.sqrt(s["n_up"]) if s["n_up"] else np.nan
        se_down = s["std_down"] / np.sqrt(s["n_down"]) if s["n_down"] else np.nan
        print(f"  SEM up = {se_up:.4f}, SEM down = {se_down:.4f}")
        print()

    if args.plot:
        import matplotlib.pyplot as plt

        pf, pm = out["pf"], out["pm"]
        fig, axes = plt.subplots(1, 2, figsize=(8, 4))
        for ax, title, means, stds, ns in (
            (
                axes[0],
                "Up state duration",
                [pf["mean_up"], pm["mean_up"]],
                [pf["std_up"], pm["std_up"]],
                [pf["n_up"], pm["n_up"]],
            ),
            (
                axes[1],
                "Down state duration",
                [pf["mean_down"], pm["mean_down"]],
                [pf["std_down"], pm["std_down"]],
                [pf["n_down"], pm["n_down"]],
            ),
        ):
            sem = [s / np.sqrt(n) if n else np.nan for s, n in zip(stds, ns)]
            ax.errorbar([1, 2], means, yerr=sem, fmt="o--")
            ax.set_xticks([1, 2])
            ax.set_xticklabels(["PF", "PM"])
            ax.set_title(title)
        fig.tight_layout()
        fig_path = args.mat.parent / "up_down_state_durations.png"
        fig.savefig(fig_path, dpi=150)
        print(f"Saved figure: {fig_path}")


if __name__ == "__main__":
    main()
