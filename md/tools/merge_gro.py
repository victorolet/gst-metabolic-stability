#!/usr/bin/env python3
"""Concatenate multiple GROMACS .gro files into one complex.gro.

Usage: merge_gro.py out.gro in1.gro [in2.gro ...]

The output atom count (line 2) is corrected automatically. Atom serial
numbers are renumbered sequentially (wrapping at 100000, per the .gro
fixed-width field). Residue numbers continue incrementing across each
input file's own numbering rather than colliding at 1. The box vector
line is carried over from the FIRST input file; downstream `gmx editconf
-bt cubic -d 1.0` redefines it anyway, so this is only a placeholder.

Order of the input files on the command line becomes the order of atoms
in the merged structure -- this must match the order of molecules you
intend to list in topol.top's [ molecules ] section exactly, since
GROMACS matches coordinates to topology purely by sequential position.
"""
import sys


def read_gro(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    title = lines[0]
    natoms = int(lines[1].strip())
    atom_lines = lines[2:2 + natoms]
    box_line = lines[2 + natoms]
    if len(atom_lines) != natoms:
        raise ValueError(f"{path}: header says {natoms} atoms, found {len(atom_lines)}")
    return title, natoms, atom_lines, box_line


def main():
    if len(sys.argv) < 3:
        print("Usage: merge_gro.py out.gro in1.gro [in2.gro ...]")
        sys.exit(1)

    out_path = sys.argv[1]
    in_paths = sys.argv[2:]

    all_atoms = []
    box_line = None
    total_atoms = 0
    resid_offset = 0

    for path in in_paths:
        title, natoms, atom_lines, box = read_gro(path)
        if box_line is None:
            box_line = box
        max_resid_this_file = 0
        for line in atom_lines:
            # .gro fixed columns: resid(5) resname(5) atomname(5) atomnum(5) x(8) y(8) z(8) [vx vy vz]
            resid = int(line[0:5])
            rest = line[5:]
            new_resid = resid + resid_offset
            max_resid_this_file = max(max_resid_this_file, resid)
            all_atoms.append(("%5d" % (new_resid % 100000)) + rest)
        resid_offset += max_resid_this_file
        total_atoms += natoms
        print(f"  {path}: {natoms} atoms, residues renumbered +{resid_offset - max_resid_this_file}")

    # .gro columns: resid(5) resname(5) atomname(5) atomnum(5) x y z ...
    # atom serial number lives at columns [15:20] -- renumber it sequentially.
    fixed_atoms = []
    for i, line in enumerate(all_atoms, 1):
        atomnum_field = "%5d" % (i % 100000)
        fixed_atoms.append(line[:15] + atomnum_field + line[20:])

    with open(out_path, "w") as fh:
        fh.write("Merged complex (protein + cofactor + ligand)\n")
        fh.write("%5d\n" % total_atoms)
        for line in fixed_atoms:
            fh.write(line + "\n")
        fh.write(box_line + "\n")

    print(f"[OK] Wrote {out_path}: {total_atoms} atoms from {len(in_paths)} input file(s)")


if __name__ == "__main__":
    main()
