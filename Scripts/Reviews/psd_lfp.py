#!/usr/bin/env python3
"""
Wavelet-based LFP power spectral density (Python port of SpectrogramAG_V4.m).

Uses the same complex-Gaussian wavelet convolution as wavelet_specgram_fast, then
averages instantaneous power |z(t)|^2 over time to obtain a 1D PSD vs frequency.

Each state is [t0 + TransientPeriod, t0 + TransientPeriod + STATE_DURATION].
STATE_DURATION is 1000 s.
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

# Starts match ARMA_*.m transients. 1000 s windows do not overlap:
# Monkey B gaps: wake|anesth 1900 s, anesth|awake 3600 s
# Monkey C gaps: wake|anesth 2000 s, anesth|awake 1250 s
STATE_DURATION = 1000.0

STATE_ORDER = ("Wakefulness", "Anesthesia", "Awakening")
STATE_TRANSIENTS = {
    "MonkeyB": {
        "Wakefulness": 100.0,
        "Anesthesia": 3000.0,
        "Awakening": 7600.0,
    },
    "MonkeyC": {
        "Wakefulness": 1000.0,
        "Anesthesia": 4000.0,
        "Awakening": 6250.0,
    },
}

# Spectrogram.m used MATLAB samples 1000–2_040_000 (~full recording, all states mixed)
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


def _as_channels_by_time(values: np.ndarray, n_time: int) -> np.ndarray:
    """Return LFP as (channels, time), handling MATLAB v7.3 transposition."""
    if values.ndim != 2:
        raise ValueError(f"MEALFP/values must be 2D, got shape {values.shape}")
    if values.shape[1] == n_time:
        return values
    if values.shape[0] == n_time:
        return values.T
    raise ValueError(
        f"MEALFP/values shape {values.shape} does not match time length {n_time}"
    )


def load_state_segments(
    mat_path: Path, monkey: str
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, float]]:
    """Load PF/PM mean LFP for each behavioral state (ARMA 250 s windows)."""
    with h5py.File(mat_path, "r") as f:
        time = np.asarray(f["MEALFP/time"]).squeeze()
        values = _as_channels_by_time(np.asarray(f["MEALFP/values"]), time.size)
        dt = float(np.asarray(f["MEALFP/dt"]).squeeze())

    t0 = float(time[0])
    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, float]] = {}
    for state in STATE_ORDER:
        transient = STATE_TRANSIENTS[monkey][state]
        t_start = t0 + transient
        t_end = t_start + STATE_DURATION
        mask = (time >= t_start) & (time <= t_end)
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            raise ValueError(
                f"{monkey} {state}: no samples in [{t_start:.1f}, {t_end:.1f}] s"
            )
        i0, i1 = int(idx[0]), int(idx[-1]) + 1
        t_vec = time[i0:i1]
        dtype = np.float32 if USE_SINGLE else np.float64
        seg_pf = np.mean(values[0:96, i0:i1], axis=0).astype(dtype)
        seg_pm = np.mean(values[96:192, i0:i1], axis=0).astype(dtype)
        out[state] = (seg_pf, seg_pm, t_vec, dt)
    return out


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


def _psd_ylim(panel_results: list[dict], f_min: float, f_max: float) -> tuple[float, float]:
    ymax = max(
        np.max(res["psd"][(res["freqs"] >= f_min) & (res["freqs"] <= f_max)])
        for res in panel_results
    )
    ymin = min(
        np.min(res["psd"][(res["freqs"] >= f_min) & (res["freqs"] <= f_max)])
        for res in panel_results
    )
    return ymin * 0.5, ymax * 3


def plot_psd(
    results: list[dict],
    f_min: float,
    f_max: float,
    output: Path | None,
) -> None:
    monkeys = list(dict.fromkeys(res["monkey"] for res in results))
    states = list(dict.fromkeys(res["state"] for res in results))
    n_row, n_col = len(states), len(monkeys)
    fig, axes = plt.subplots(
        n_row, n_col, figsize=(5 * n_col, 3.8 * n_row), squeeze=False
    )

    for col, monkey in enumerate(monkeys):
        monkey_results = [res for res in results if res["monkey"] == monkey]
        ylim = _psd_ylim(monkey_results, f_min, f_max)

        for row, state in enumerate(states):
            ax = axes[row, col]
            panel = [res for res in monkey_results if res["state"] == state]
            for res in panel:
                mask = (res["freqs"] >= f_min) & (res["freqs"] <= f_max)
                ax.loglog(
                    res["freqs"][mask],
                    res["psd"][mask],
                    lw=2.0,
                    color=AREA_COLORS[res["area"]],
                    label=AREA_LABELS[res["area"]],
                )

            ax.set_xlim(f_min, f_max)
            ax.set_ylim(*ylim)
            ax.grid(True, which="both", alpha=0.25)
            if row == 0:
                ax.set_title(MONKEY_LABELS[monkey])
            if col == 0:
                ax.set_ylabel(f"{state}\nPower spectral density (a.u.)")
            if row == n_row - 1:
                ax.set_xlabel("Frequency (Hz)")
            if row == 0:
                ax.legend(loc="upper right", framealpha=0.9)

    fig.suptitle("LFP power spectral density", y=1.01, fontsize=12)
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
        default=Path(__file__).resolve().parent / "figures" / "lfp_psd_loglog_by_state.png",
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
        state_segments = load_state_segments(mat_path, name)

        for state in STATE_ORDER:
            seg_pf, seg_pm, t_vec, dt = state_segments[state]
            segments = {"PF": seg_pf, "PM": seg_pm}
            transient = STATE_TRANSIENTS[name][state]
            print(
                f"  {state}: t0+{transient:g}–{transient + STATE_DURATION:g} s "
                f"({t_vec[0]:.1f}–{t_vec[-1]:.1f} s absolute, "
                f"{seg_pf.size} samples, dt={dt:g} s, fs={1 / dt:.1f} Hz)"
            )

            for area_key, _area_label in areas:
                label = f"{name} {state} {area_key}"
                print(f"    computing wavelet PSD: {label}")
                psd = wavelet_psd(segments[area_key], freqs, dt)
                results.append(
                    {
                        "label": label,
                        "monkey": name,
                        "state": state,
                        "area": area_key,
                        "freqs": freqs,
                        "psd": psd,
                    }
                )

                slope = fit_loglog_slope(freqs, psd, 2.0, 40.0)
                print(f"    1/f slope (2–40 Hz): {slope:.2f}")
                for band, (f_lo, f_hi) in BANDS.items():
                    if f_hi <= args.f_max:
                        print(
                            f"    {band:5s} power ({f_lo:g}–{f_hi:g} Hz): "
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
