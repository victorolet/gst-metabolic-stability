#!/usr/bin/env python3
"""Merge per-frame DELTA TOTAL binding-energy values from multiple
gmx_MMPBSA chunk runs (produced by run_mmpbsa_array.sh) into one pooled
mean/SEM, instead of averaging each chunk's own mean.

Reads every chunk_*/FINAL_RESULTS_MMPBSA.csv under the given directory,
locates the "DELTA Energy Terms" section, pulls the (frame, TOTAL) pairs,
concatenates across chunks (flagging any frame that appears in more than
one chunk -- should not happen if chunk ranges were computed correctly),
and reports the pooled binding energy estimate.

Usage: merge_mmpbsa_chunks.py <chunks-dir> [--out merged_delta_total.csv]
"""
import argparse
import csv
import glob
import os
import statistics
import sys


def parse_delta_total(csv_path):
    """Return list of (frame, total) from the DELTA Energy Terms section
    of a gmx_MMPBSA -eo CSV file."""
    with open(csv_path, newline="") as fh:
        rows = list(csv.reader(fh))

    section_start = None
    for i, row in enumerate(rows):
        joined = ",".join(row).upper()
        if "DELTA" in joined and "ENERGY" in joined and "TERM" in joined:
            section_start = i
            break
    if section_start is None:
        raise ValueError(f"no 'DELTA Energy Terms' section header found in {csv_path}")

    # header row is the next non-empty row after the section title
    header_idx = None
    for i in range(section_start + 1, len(rows)):
        if any(cell.strip() for cell in rows[i]):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"no header row found after DELTA section in {csv_path}")

    header = [c.strip().upper() for c in rows[header_idx]]
    frame_col = next((j for j, c in enumerate(header) if "FRAME" in c), None)
    total_col = next((j for j, c in enumerate(header) if "TOTAL" in c), None)
    if frame_col is None or total_col is None:
        raise ValueError(
            f"could not find FRAME/TOTAL columns in {csv_path} header: {header}"
        )

    out = []
    for row in rows[header_idx + 1:]:
        if not any(cell.strip() for cell in row):
            break  # blank line = end of section
        try:
            frame = row[frame_col].strip()
            total = float(row[total_col])
        except (ValueError, IndexError):
            break
        out.append((frame, total))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("chunks_dir", help="directory containing chunk_*/FINAL_RESULTS_MMPBSA.csv")
    ap.add_argument("--out", default=None, help="write merged per-frame CSV here")
    args = ap.parse_args()

    csv_paths = sorted(glob.glob(os.path.join(args.chunks_dir, "chunk_*", "FINAL_RESULTS_MMPBSA.csv")))
    if not csv_paths:
        sys.exit(f"ERROR: no chunk_*/FINAL_RESULTS_MMPBSA.csv files found under {args.chunks_dir}")

    all_rows = []      # (chunk_label, frame, total)
    seen_frames = {}   # frame -> chunk_label, to flag accidental overlap
    for p in csv_paths:
        chunk_label = os.path.basename(os.path.dirname(p))
        try:
            pairs = parse_delta_total(p)
        except ValueError as e:
            print(f"WARNING: skipping {p}: {e}", file=sys.stderr)
            continue
        if not pairs:
            print(f"WARNING: {p} parsed but yielded 0 frames", file=sys.stderr)
        for frame, total in pairs:
            if frame in seen_frames:
                print(
                    f"WARNING: frame {frame} appears in both {seen_frames[frame]} and "
                    f"{chunk_label} -- check for overlapping chunk ranges",
                    file=sys.stderr,
                )
            seen_frames[frame] = chunk_label
            all_rows.append((chunk_label, frame, total))

    if not all_rows:
        sys.exit("ERROR: parsed 0 frames total across all chunks -- nothing to merge")

    totals = [t for _, _, t in all_rows]
    n = len(totals)
    mean = statistics.mean(totals)
    stdev = statistics.stdev(totals) if n > 1 else 0.0
    sem = stdev / (n ** 0.5) if n > 1 else 0.0

    print(f"Chunks merged: {len(csv_paths)}")
    print(f"Total frames:  {n}")
    print(f"DELTA TOTAL (pooled): {mean:.4f} +/- {sem:.4f} kcal/mol (SEM), stdev {stdev:.4f}")

    if args.out:
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["chunk", "frame", "DELTA_TOTAL"])
            w.writerows(all_rows)
        print(f"[OK] wrote per-frame merged table to {args.out}")


if __name__ == "__main__":
    main()
