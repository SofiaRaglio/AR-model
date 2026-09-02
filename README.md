# AR-model

Analysis code for the paper *"Local field potentials predictability from spiking activity as a probe of network complexity"*.

The repository trains autoregressive (AR) models to predict local field potentials (LFP) from multi-unit activity (MUA) in prefrontal (PF) and premotor (PM) cortex, across wakefulness, anaesthesia, and awakening, in two macaques (Monkey B and Monkey C).

## Repository layout

```
Scripts/
├── *.m                 # Common scripts (shared across both monkeys)
├── MonkeyB/
│   ├── *.m             # Monkey B–specific analyses
│   └── Kernel/         # Precomputed AR kernels and performance metrics
├── MonkeyC/
│   ├── *.m             # Monkey C–specific analyses
│   └── Kernel/         # Precomputed AR kernels and performance metrics
└── Reviews/            # Additional analyses requested by reviewers (Python)
```

Scripts in `Scripts/` (outside the monkey folders) are **common** to both subjects. Each monkey folder contains the same analysis pipeline, adapted to subject-specific channels, transient periods, and file paths.

Original electrophysiology recordings (`MEAMUALFP.mat`, `PreProcData.mat`) are **not** included. Same for their fragmented version used for the percolation analysis. Update the `path/to/...` placeholders in the MATLAB scripts before running.

## Behavioral states

Each state is defined by a transient period at the start of the training window:

| State        | Monkey B | Monkey C |
|--------------|----------|----------|
| Wakefulness  | 100 s    | 1000 s   |
| Anaesthesia  | 3000 s   | 4000 s   |
| Awakening    | 7600 s   | 6250 s   |

Precomputed results live in `MonkeyB/Kernel/` and `MonkeyC/Kernel/`. Files are named by area and transient, e.g.:

- `KernelAndPerfPF_trans100.mat` — intra-area PF
- `KernelAndPerfPM_trans3000.mat` — intra-area PM
- `KernelAndPerfPF2PM_trans7600.mat` — PF → PM (cross-area)
- `KernelAndPerfPM2PF_trans100.mat` — PM → PF (cross-area)

Run the `ARMA_*.m` scripts to retrain kernels from raw data; use the plotting and analysis scripts to reproduce paper figures from the saved `.mat` files.

## Common scripts (`Scripts/`)

| Script | Paper figure | Description |
|--------|--------------|-------------|
| `Spectrogram.m` | Fig. 1a, Suppl. Fig. 1a | LFP spectrogram for PF and PM |
| `detectingUpAndDownStates.m` | Fig. 1c | Up/down state detection |
| `Boxplot.m` | Fig. 2b | Intra-area kernel performance across states |
| `BoxplotInterArea.m` | Fig. 5a | Inter-area kernel performance across states |
| `Test_Percolation_PrePost.m` | Fig. 5c | Percolation around state transitions |

## Monkey-specific scripts (`Scripts/MonkeyB/`, `Scripts/MonkeyC/`)

| Script | Paper figure | Description |
|--------|--------------|-------------|
| `ARMA_PF.m`, `ARMA_PM.m` | Fig. 2b | Train intra-area AR models |
| `ARMA_PFtoPM.m`, `ARMA_PMtoPF.m` | Fig. 5a | Train cross-area AR models |
| `plotStarsInPCn.m` | Fig. 3a,c,d (B); Suppl. Fig. 2c (C) | Kernel PCA visualisation |
| `KernelDistDependencePF.m`, `KernelDistDependencePM.m` | Fig. 3d | Kernel distance dependence of contribution to reconstruction |
| `PcaKernels.m` | Fig. 4a (B); Suppl. Fig. 3a (C) | Correlations between kernel PCA components |
| `PC_trunc.m` | Fig. 4b (B); Fig. 2b (C) | PCA truncation at 90% variance |
| `Percolation.m` | — | Percolation function |
| `Run_Percolation.m` / `Run_Percolation_C.m` | — | Percolation on correlation graphs |
| `TransitionsPrePost.m` | Fig. 6 (B); Suppl. Fig. 5 (C) | Up-to-down transition analysis |

`Run_Percolation.m` expects per-condition `.mat` files (`Wakefulness.mat`, `Anaesthesia.mat`, `Awakening.mat`) in the working directory.

## Reviewer analyses (`Scripts/Reviews/`)

Python scripts with supplementary analyses requested during peer review. Several are ports of the common MATLAB figure scripts.

| Script | Paper Figure | Description |
|--------|--------------|-------------|
| `boxplot_kernels.py` | | Port of `Boxplot.m` |
| `boxplot_interarea_kernels.py` | | Port of `BoxplotInterArea.m` |
| `detecting_up_and_down_states.py` | | Port of `detectingUpAndDownStates.m` |
| `psd_lfp.py` | Suppl. Fig. 1b | Wavelet-based LFP power spectral density |
| `percolation_analysis.py` | | Shared percolation analysis module |
| `run_percolation_summary.py` | | Run percolation for both monkeys; export summary table |
| `evaluate_kernel_metrics.py` | | Evaluate precomputed kernels on held-out test windows |
| `arma_average_lfp.py` | Suppl. Fig. 2 | Compute the performance comparison between average to average prediction, per-channel prediction and average subtracted per-channel prediction |


## Requirements

- **MATLAB** — common and monkey-specific `.m` scripts
- **Python 3** with `numpy`, `scipy`, `matplotlib`, and `h5py` — `Scripts/Reviews/`
