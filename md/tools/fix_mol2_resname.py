#!/usr/bin/env python3
"""Force every atom's substructure (residue) name in a .mol2 file to a
fixed value, in place. Replaces Harry's manual step ("Edit ligand_fixed.mol2
removing '*****' and replacing with the ligand name, i.e. 'UNL'") -- rather
than hunting for a specific placeholder string, this rewrites the actual
subst_name field of every @<TRIPOS>ATOM line and the name in
@<TRIPOS>SUBSTRUCTURE, whatever the upstream tool (obabel, meeko, CGenFF's
own export, ...) happened to put there.

Usage: fix_mol2_resname.py in.mol2 out.mol2 RESNAME
"""
import sys


def fix_atom_line(line, resname):
    # TRIPOS ATOM record: id name x y z type [subst_id [subst_name [charge]]]
    parts = line.split()
    if len(parts) < 6:
        return line  # malformed/blank -- leave untouched
    while len(parts) < 8:
        parts.append({6: "1", 7: resname}.get(len(parts), "0"))
    parts[7] = resname
    # rebuild with mol2's conventional column widths (7 int, then fields)
    return "%7s %-8s %9s %9s %9s %-6s %4s %-8s %9s\n" % (
        parts[0], parts[1], parts[2], parts[3], parts[4],
        parts[5], parts[6], parts[7], parts[8] if len(parts) > 8 else "0.0000",
    )


def main():
    if len(sys.argv) != 4:
        print("Usage: fix_mol2_resname.py in.mol2 out.mol2 RESNAME")
        sys.exit(1)
    in_path, out_path, resname = sys.argv[1], sys.argv[2], sys.argv[3]

    with open(in_path) as fh:
        lines = fh.readlines()

    out = []
    section = None
    n_fixed = 0
    for line in lines:
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            out.append(line)
            continue
        if section == "@<TRIPOS>ATOM" and line.strip():
            out.append(fix_atom_line(line, resname))
            n_fixed += 1
        elif section == "@<TRIPOS>SUBSTRUCTURE" and line.strip():
            parts = line.split()
            if len(parts) >= 2:
                parts[1] = resname
            out.append(" ".join(parts) + "\n")
        else:
            out.append(line)

    if n_fixed == 0:
        raise SystemExit(f"ERROR: no @<TRIPOS>ATOM records found/fixed in {in_path} "
                          "-- is this a valid mol2 file?")

    with open(out_path, "w") as fh:
        fh.writelines(out)
    print(f"[OK] Renamed residue on {n_fixed} atoms -> '{resname}' in {out_path}")


if __name__ == "__main__":
    main()
