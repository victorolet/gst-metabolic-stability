#!/usr/bin/env python3
"""Render a quick PNG line plot from a GROMACS-style .xvg file (columns:
x, y[, y2, ...]). Reads @ title / @ xaxis label / @ yaxis label / @ sN
legend lines when present for automatic labelling; CLI flags override.

Usage: plot_xvg.py input.xvg output.png [--title T] [--xlabel X] [--ylabel Y]
"""
import argparse
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_xvg(path):
    title = xlabel = ylabel = None
    legends = {}
    data = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("@"):
                m = re.search(r'title\s+"(.*)"', line)
                if m:
                    title = m.group(1)
                m = re.search(r'xaxis\s+label\s+"(.*)"', line)
                if m:
                    xlabel = m.group(1)
                m = re.search(r'yaxis\s+label\s+"(.*)"', line)
                if m:
                    ylabel = m.group(1)
                m = re.search(r's(\d+)\s+legend\s+"(.*)"', line)
                if m:
                    legends[int(m.group(1))] = m.group(2)
                continue
            if line.startswith("#") or not line.strip():
                continue
            data.append([float(v) for v in line.split()])
    return title, xlabel, ylabel, legends, data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xvg")
    ap.add_argument("png")
    ap.add_argument("--title", default=None)
    ap.add_argument("--xlabel", default=None)
    ap.add_argument("--ylabel", default=None)
    args = ap.parse_args()

    title, xlabel, ylabel, legends, data = read_xvg(args.xvg)
    if not data:
        sys.exit(f"ERROR: no data rows parsed from {args.xvg}")

    ncols = len(data[0])
    xs = [row[0] for row in data]

    plt.figure(figsize=(7, 4.5))
    for col in range(1, ncols):
        ys = [row[col] for row in data]
        label = legends.get(col - 1, f"col{col}") if ncols > 2 else None
        plt.plot(xs, ys, label=label, linewidth=1.2)

    plt.title(args.title or title or args.xvg)
    plt.xlabel(args.xlabel or xlabel or "")
    plt.ylabel(args.ylabel or ylabel or "")
    if ncols > 2:
        plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.png, dpi=150)
    print(f"[OK] wrote {args.png}")


if __name__ == "__main__":
    main()
