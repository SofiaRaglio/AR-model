#!/usr/bin/env python3
"""Run percolation analysis for Monkey B and Monkey C; print combined summary table."""

from pathlib import Path

from percolation_analysis import (
    MONKEY_B_CONFIG,
    MONKEY_C_CONFIG,
    print_summary_table,
    run_subject,
    summary_table_rows,
)

SUMMARY_DIR = Path(__file__).resolve().parent


def main():
    all_rows = []
    results = {}

    for cfg in (MONKEY_B_CONFIG, MONKEY_C_CONFIG):
        cfg.save_figures = True
        out = run_subject(cfg)
        results[cfg.name] = out
        all_rows.extend(summary_table_rows(cfg.name, out["stats"]))

    print_summary_table(all_rows)

    csv_path = SUMMARY_DIR / "percolation_summary_both_animals.csv"
    import csv

    fieldnames = [
        "subject",
        "comparison",
        "mean_diff",
        "wilcoxon_p",
        "wilcoxon_p_bonf",
        "ttest_p_bonf",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nCSV saved to {csv_path}")

    return results, all_rows


if __name__ == "__main__":
    main()
