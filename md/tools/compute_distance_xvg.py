#!/usr/bin/env python3
"""Read two GROMACS .xvg files (each: time, x, y, z -- the output of
`gmx traj -com -ox` on a single-atom index group) and write a third .xvg
with the per-frame Euclidean distance between them, in nm (GROMACS's native
length unit).

Usage: compute_distance_xvg.py atom_a.xvg atom_b.xvg out_distance.xvg
"""
import math
import sys


def read_xvg_xyz(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("@", "#")):
                continue
            parts = line.split()
            t, x, y, z = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
            rows.append((t, x, y, z))
    return rows


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: compute_distance_xvg.py atom_a.xvg atom_b.xvg out_distance.xvg")
    a_path, b_path, out_path = sys.argv[1:4]

    a = read_xvg_xyz(a_path)
    b = read_xvg_xyz(b_path)
    if not a or not b:
        sys.exit(f"ERROR: no data rows parsed from {a_path if not a else b_path}")
    if len(a) != len(b):
        sys.exit(f"ERROR: frame count mismatch ({len(a)} vs {len(b)}) between {a_path} and {b_path}")

    with open(out_path, "w") as fh:
        fh.write('@    title "Ligand warhead - GSH thiolate (SG2) distance"\n')
        fh.write('@    xaxis label "Time (ns)"\n')
        fh.write('@    yaxis label "Distance (nm)"\n')
        fh.write("@TYPE xy\n")
        for (ta, xa, ya, za), (tb, xb, yb, zb) in zip(a, b):
            if abs(ta - tb) > 1e-6:
                sys.exit(f"ERROR: time mismatch between inputs at t={ta} vs t={tb} -- frames not aligned")
            d = math.sqrt((xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2)
            fh.write(f"{ta:.4f}  {d:.4f}\n")

    print(f"[OK] wrote {out_path} ({len(a)} frames)")


if __name__ == "__main__":
    main()
