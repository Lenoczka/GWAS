#!/usr/bin/env python3

import sys
import os
import math

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def normalize_chr(value):
    s = str(value).strip().upper()
    if s.startswith("CHR"):
        s = s[3:]
    mapping = {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}
    if s in mapping:
        return mapping[s]
    try:
        return int(float(s))
    except ValueError:
        return np.nan


def read_assoc_table(path):
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df = df.loc[[x == 'ADD' for x in df['TEST']]].copy()

    required = {"CHR", "BP", "P"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Input file must contain columns {sorted(required)}; missing: {sorted(missing)}"
        )

    df = df[["CHR", "BP", "P"]].copy()
    df["CHR"] = df["CHR"].apply(normalize_chr)
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
    df["P"] = pd.to_numeric(df["P"], errors="coerce")

    df = df.dropna(subset=["CHR", "BP", "P"])
    df = df[(df["P"] > 0) & (df["P"] <= 1)]
    df["CHR"] = df["CHR"].astype(int)

    if df.empty:
        raise ValueError("No valid rows remained after filtering invalid CHR/BP/P values.")

    return df


def make_qq_plot(pvals, out_png):
    pvals = np.sort(np.asarray(pvals))
    n = len(pvals)

    observed = -np.log10(pvals)
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)

    maxval_exp = expected.max() * 1.05
    maxval_obs = observed.max() * 1.05

    plt.figure(figsize=(6, 6))
    plt.scatter(expected, observed, s=8, alpha=0.7)
    plt.plot([0, maxval_exp], [0, maxval_exp], linestyle="--", linewidth=1)
    plt.xlabel("Expected -log10(P)")
    plt.ylabel("Observed -log10(P)")
    plt.title("Q-Q plot")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def make_manhattan_plot(df, out_png):
    df = df.sort_values(["CHR", "BP"]).copy()
    df["minus_log10_p"] = -np.log10(df["P"])

    chromosomes = sorted(df["CHR"].unique())
    x_positions = []
    tick_positions = []
    tick_labels = []
    offset = 0

    for i, chrom in enumerate(chromosomes):
        chr_df = df[df["CHR"] == chrom].sort_values("BP").copy()
        chr_df["x"] = chr_df["BP"] + offset
        x_positions.append(chr_df)

        tick_positions.append((chr_df["x"].min() + chr_df["x"].max()) / 2)
        tick_labels.append(str(chrom) if chrom <= 22 else {23: "X", 24: "Y", 25: "XY", 26: "MT"}.get(chrom, str(chrom)))

        offset = chr_df["x"].max()

    plot_df = pd.concat(x_positions, ignore_index=True)

    plt.figure(figsize=(12, 5))
    for i, chrom in enumerate(chromosomes):
        chr_df = plot_df[plot_df["CHR"] == chrom]
        plt.scatter(
            chr_df["x"],
            chr_df["minus_log10_p"],
            s=8,
            alpha=0.75,
            label=str(chrom) if i < 2 else None,
        )

    if (df["P"] < 5e-8).any():
        plt.axhline(-math.log10(5e-8), linestyle="--", linewidth=1)

    plt.xticks(tick_positions, tick_labels, rotation=0)
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(P)")
    plt.title("Manhattan plot")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    if len(sys.argv) != 2:
        sys.stderr.write(
            f"Usage: {os.path.basename(sys.argv[0])} <plink_assoc_results.txt>\n"
        )
        sys.exit(1)

    in_path = sys.argv[1]
    base = os.path.splitext(os.path.basename(in_path))[0]

    qq_png = f"{base}.qq.png"
    manhattan_png = f"{base}.manhattan.png"

    df = read_assoc_table(in_path)
    make_qq_plot(df["P"].values, qq_png)
    make_manhattan_plot(df, manhattan_png)

    print(f"Wrote: {qq_png}")
    print(f"Wrote: {manhattan_png}")


if __name__ == "__main__":
    main()
