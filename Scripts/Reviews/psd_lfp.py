#!/usr/bin/env python3
"""
Wavelet-based LFP power spectral density (Python port of SpectrogramAG_V4.m).

Uses the same complex-Gaussian wavelet convolution as wavelet_specgram_fast, then
averages instantaneous power |z(t)|^2 over time to obtain a 1D PSD vs frequency.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

SCRIPTS_ROOT = Path(__file__).resolve().parent.parent

DEFAULT_DATASETS = {
    "MonkeyB": SCRIPTS_ROOT / "MonkeyB" / "MEAMUALFP.mat",
    "MonkeyC": SCRIPTS_ROOT / "MonkeyC" / "MEAMUALFP.mat",
}

MONKEY_LABELS = {
    "MonkeyB": "Monkey B",
    "MonkeyC": "Monkey C",
}

AREA_COLORS = {
    "PF": "#00BCD4",  # cyan — prefrontal
    "PM": "#1565C0",  # blue — premotor
}

AREA_LABELS = {
    "PF": "Prefrontal",
    "PM": "Premotor",
}

# SpectrogramAG_V4.m defaults
DT = 0.005
TIME_START_IDX = 1000       # MATLAB 1-based index
TIME_END_IDX = 2_040_000    # MATLAB 1-based index
BANDWIDTH = 0.05
SD_TIMES = 4.0
USE_SINGLE = True

BANDS = {
    "delta": (1.0, 4.0),
    "theta": (4.0, 8.0),
    "alpha": (8.0, 12.0),
    "beta": (13.0, 30.0),
    "gamma": (30.0, 80.0),
}


def load_segment_ag_v4(mat_path: Path) -> tuple[np.ndarray, np.ndarray, float]:
    """Load PF/PM mean LFP traces using SpectrogramAG_V4 indexing and channel groups."""
    with h5py.File(mat_path, "r") as f:
        time = np.asarray(f["MEALFP/time"]).squeeze()
        values = np.asarray(f["MEALFP/values"]).T  # channels x time
        dt = float(np.asarray(f["MEALFP/dt"]).squeeze())

    t0 = time[0]
    i_start = max(0, TIME_START_IDX - 1)
    i_end = min(values.shape[1], TIME_END_IDX)
    t_start = time[i_start]
    t_end = time[i_end - 1]

    time_s = max(0, int(round((t_start - t0) / dt)))
    time_e = min(values.shape[1], int(round((t_end - t0) / dt)) + 1)

    seg_pf = np.mean(values[0:96, time_s:time_e], axis=0).astype(np.float64)
    seg_pm = np.mean(values[96:192, time_s:time_e], axis=0).astype(np.float64)
    t_vec = time[time_s:time_e]

    if USE_SINGLE:
        seg_pf = seg_pf.astype(np.float32)
        seg_pm = seg_pm.astype(np.float32)

    return seg_pf, seg_pm, t_vec, dt


def frequency_centers_log(f1: float, f2: float, n_freq: int) -> np.ndarray:
    """Log-spaced frequency bin centers (geometric mean of log-spaced edges)."""
    f_edges = np.exp(np.linspace(np.log(f1), np.log(f2), n_freq + 1))
    return np.sqrt(f_edges[:-1] * f_edges[1:])


def wavelet_psd(
    x: np.ndarray,
    freqs: np.ndarray,
    dt: float,
    rel_bandwidth: float = BANDWIDTH,
    sd_times: float = SD_TIMES,
    use_single: bool = USE_SINGLE,
) -> np.ndarray:
    """
    Wavelet PSD: mean instantaneous power |conv(x, WL)|^2 over time.

    Port of the per-frequency loop in wavelet_specgram_fast (SpectrogramAG_V4.m),
    without the moving-average temporal downsampling.
    """
    x = np.asarray(x, dtype=np.float32 if use_single else np.float64).ravel()
    n = x.size
    n_freq = freqs.size
    psd = np.empty(n_freq, dtype=np.float64)

    sigma_min = 1.0 / (2.0 * np.pi * freqs.min() * rel_bandwidth)
    max_half_len = int(np.ceil(sd_times * sigma_min / dt))
    max_wl = 2 * max_half_len + 1
    n_fft = int(2 ** np.ceil(np.log2(n + max_wl - 1)))

    x_fft = np.fft.fft(x, n_fft)

    for k, fk in enumerate(freqs):
        sigma = 1.0 / (2.0 * np.pi * fk * rel_bandwidth)
        half_len = int(np.ceil(sd_times * sigma / dt))
        xx = np.arange(-half_len, half_len + 1, dtype=np.float64) * dt

        wl = (
            (1.0 / np.sqrt(sigma * np.sqrt(np.pi)))
            * np.exp(-(xx**2) / (2.0 * sigma**2))
            * np.exp(1j * 2.0 * np.pi * fk * xx)
        )
        if use_single:
            wl = wl.astype(np.complex64)

        y_fft = np.fft.fft(wl, n_fft)
        z = np.fft.ifft(x_fft * y_fft) * dt
        p = np.abs(z[half_len : half_len + n]) ** 2
        psd[k] = float(np.mean(p))

        if (k + 1) % 20 == 0 or k + 1 == n_freq:
            print(f"    frequency {k + 1}/{n_freq} ({fk:.3g} Hz)")

    return psd


def band_power(freqs: np.ndarray, psd: np.ndarray, f_lo: float, f_hi: float) -> float:
    mask = (freqs >= f_lo) & (freqs <= f_hi)
    return float(np.trapezoid(psd[mask], freqs[mask]))


def fit_loglog_slope(
    freqs: np.ndarray, psd: np.ndarray, f_lo: float, f_hi: float
) -> float:
    mask = (freqs >= f_lo) & (freqs <= f_hi) & (psd > 0)
    slope, _ = np.polyfit(np.log10(freqs[mask]), np.log10(psd[mask]), 1)
    return float(slope)


def plot_psd(
    results: list[dict],
    f_min: float,
    f_max: float,
    output: Path | None,
) -> None:
    monkeys = list(dict.fromkeys(res["monkey"] for res in results))
    fig, axes = plt.subplots(1, len(monkeys), figsize=(5 * len(monkeys), 4.5), squeeze=False)

    for col, monkey in enumerate(monkeys):
        ax = axes[0, col]
        monkey_results = [res for res in results if res["monkey"] == monkey]

        for res in monkey_results:
            mask = (res["freqs"] >= f_min) & (res["freqs"] <= f_max)
            ax.loglog(
                res["freqs"][mask],
                res["psd"][mask],
                lw=2.0,
                color=AREA_COLORS[res["area"]],
                label=AREA_LABELS[res["area"]],
            )

        ymax = max(
            np.max(res["psd"][(res["freqs"] >= f_min) & (res["freqs"] <= f_max)])
            for res in monkey_results
        )
        ymin = min(
            np.min(res["psd"][(res["freqs"] >= f_min) & (res["freqs"] <= f_max)])
            for res in monkey_results
        )

        ax.set_xlim(f_min, f_max)
        ax.set_ylim(ymin * 0.5, ymax * 3)
        ax.set_title(MONKEY_LABELS[monkey])
        ax.legend(loc="upper right", framealpha=0.9)
        ax.grid(True, which="both", alpha=0.25)
        if col == 0:
            ax.set_ylabel("Power spectral density (a.u.)")
        ax.set_xlabel("Frequency (Hz)")

    fig.suptitle("LFP power spectral density", y=1.02, fontsize=12)
    fig.tight_layout()

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=200, bbox_inches="tight")
        print(f"Saved figure to {output}")
        pdf_path = output.with_suffix(".pdf")
        fig.savefig(pdf_path, bbox_inches="tight")
        print(f"Saved figure to {pdf_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=list(DEFAULT_DATASETS),
        choices=list(DEFAULT_DATASETS),
    )
    parser.add_argument(
        "--area",
        choices=("PF", "PM", "both"),
        default="both",
        help="Prefrontal (ch 1:96) and/or premotor (ch 97:192) [both]",
    )
    parser.add_argument(
        "--f-min",
        type=float,
        default=0.5,
        help="Minimum frequency (Hz) [0.5]",
    )
    parser.add_argument(
        "--f-max",
        type=float,
        default=100.0,
        help="Maximum frequency (Hz); Nyquist at 200 Hz sampling is 100 [100]",
    )
    parser.add_argument(
        "--n-freq",
        type=int,
        default=80,
        help="Number of log-spaced frequency samples [80]",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "figures" / "lfp_psd_loglog.png",
    )
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    freqs = frequency_centers_log(args.f_min, args.f_max, args.n_freq)
    areas: list[tuple[str, str]] = []
    if args.area in ("PF", "both"):
        areas.append(("PF", "Prefrontal (ch 1:96)"))
    if args.area in ("PM", "both"):
        areas.append(("PM", "Premotor (ch 97:192)"))

    results: list[dict] = []
    for name in args.datasets:
        mat_path = DEFAULT_DATASETS[name]
        print(f"\n=== {MONKEY_LABELS[name]} ===")
        seg_pf, seg_pm, t_vec, dt = load_segment_ag_v4(mat_path)
        segments = {"PF": seg_pf, "PM": seg_pm}

        print(
            f"  segment: {t_vec[0]:.1f}–{t_vec[-1]:.1f} s "
            f"({seg_pf.size} samples, dt={dt:g} s, fs={1/dt:.1f} Hz)"
        )

        for area_key, area_label in areas:
            label = f"{name} {area_key}"
            print(f"  computing wavelet PSD: {label}")
            psd = wavelet_psd(segments[area_key], freqs, dt)
            results.append(
                {
                    "label": label,
                    "monkey": name,
                    "area": area_key,
                    "freqs": freqs,
                    "psd": psd,
                }
            )

            slope = fit_loglog_slope(freqs, psd, 2.0, 40.0)
            print(f"  1/f slope (2–40 Hz): {slope:.2f}")
            for band, (f_lo, f_hi) in BANDS.items():
                if f_hi <= args.f_max:
                    print(
                        f"  {band:5s} power ({f_lo:g}–{f_hi:g} Hz): "
                        f"{band_power(freqs, psd, f_lo, f_hi):.4g}"
                    )

    if not args.no_show:
        plot_psd(results, f_min=args.f_min, f_max=args.f_max, output=args.output)
        plt.show()
    else:
        plt.ioff()
        plot_psd(results, f_min=args.f_min, f_max=args.f_max, output=args.output)
        plt.close("all")


if __name__ == "__main__":
    main()
