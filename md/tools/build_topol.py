#!/usr/bin/env python3
"""Assemble a complex topol.top from a protein-only topol.top (as produced
by `gmx pdb2gmx`) plus one or more extra molecules (GSH, ligand), replacing
the manual "add unl.prm, unl.itp, GSH.prm, and GSH.itp to topol.top" step
in Harry's workflow.

What this does, precisely:
  1. Inserts `#include "<name>.prm"` lines for every extra molecule right
     after the `forcefield.itp` include (parameters must be visible before
     any moleculetype that uses them).
  2. Inserts `#include "<name>.itp"` lines (each followed by an `#ifdef
     POSRES` block referencing `posre_<name>.itp`) right before the water
     topology include -- this matches the order GROMACS requires (all
     non-water molecule types before water) and mirrors Harry's own
     topol.top example.
  3. Appends `<MOLNAME>   1` lines to `[ molecules ]`, in the same order
     the molecules will appear in the merged .gro file (this must match
     merge_gro.py's input order exactly).

The posre_<name>.itp files referenced inside the #ifdef POSRES blocks do
not need to exist yet when this script runs -- ions.mdp/minim.mdp don't
define POSRES, so grompp skips those blocks for EM; by the time NVT runs
(which does define POSRES via nvt.mdp's `define = -DPOSRES`), genrestr
will have already produced them.

Usage:
  build_topol.py protein_topol.top out_topol.top \\
      --molecule GSH:GSH.itp:GSH.prm:posre_GSH.itp \\
      --molecule UNL:unl.itp:unl.prm:posre_unl.itp
"""
import argparse
import re


def parse_molecule_arg(s):
    parts = s.split(":")
    if len(parts) != 4:
        raise argparse.ArgumentTypeError(
            f"--molecule must be NAME:ITP:PRM:POSRE_ITP, got: {s}"
        )
    name, itp, prm, posre = parts
    return {"name": name, "itp": itp, "prm": prm, "posre": posre}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("protein_topol")
    ap.add_argument("out_topol")
    ap.add_argument("--molecule", action="append", type=parse_molecule_arg, required=True,
                     help="NAME:ITP:PRM:POSRE_ITP -- may be repeated, in .gro order")
    args = ap.parse_args()

    with open(args.protein_topol) as fh:
        lines = fh.readlines()

    forcefield_re = re.compile(r'#include\s+".*forcefield\.itp"')
    water_re = re.compile(r'#include\s+".*(spc|spce|tip3p|tip4p)\.itp"', re.IGNORECASE)
    molecules_header_re = re.compile(r'^\s*\[\s*molecules\s*\]')

    ff_idx = next((i for i, l in enumerate(lines) if forcefield_re.search(l)), None)
    water_idx = next((i for i, l in enumerate(lines) if water_re.search(l)), None)
    mol_header_idx = next((i for i, l in enumerate(lines) if molecules_header_re.match(l)), None)

    missing = [name for name, idx in
               [("forcefield.itp include", ff_idx),
                ("water topology include", water_idx),
                ("[ molecules ] section", mol_header_idx)]
               if idx is None]
    if missing:
        raise SystemExit(
            "ERROR: could not find expected anchor(s) in "
            f"{args.protein_topol}: {', '.join(missing)}. "
            "This script expects a topol.top as produced by `gmx pdb2gmx` "
            "-- if the format looks different, the anchors below may need "
            "adjusting rather than trusting this output blindly."
        )
    if not (ff_idx < water_idx < mol_header_idx):
        raise SystemExit(
            "ERROR: anchor lines found in an unexpected order "
            f"(forcefield.itp@{ff_idx}, water@{water_idx}, [molecules]@{mol_header_idx}) "
            "-- refusing to guess where to insert content."
        )

    # --- Block 1: parameter includes, right after forcefield.itp ---
    prm_block = ["\n; extra molecule parameters (added by build_topol.py)\n"]
    for mol in args.molecule:
        prm_block.append(f'#include "{mol["prm"]}"\n')

    # --- Block 2: molecule itp + posre, right before the water include ---
    itp_block = ["\n; extra molecule topologies (added by build_topol.py)\n"]
    for mol in args.molecule:
        itp_block.append(f'#include "{mol["itp"]}"\n')
        itp_block.append("#ifdef POSRES\n")
        itp_block.append(f'#include "{mol["posre"]}"\n')
        itp_block.append("#endif\n")
        itp_block.append("\n")

    # --- Block 3: [ molecules ] entries, appended at the end of the file ---
    mol_lines = ["; extra molecules (added by build_topol.py)\n"]
    for mol in args.molecule:
        mol_lines.append(f'{mol["name"]:<16s} 1\n')

    # Assemble, working from the bottom up so earlier indices stay valid.
    out = list(lines)
    out = out[:water_idx] + itp_block + out[water_idx:]
    # water_idx shifted by len(itp_block) items inserted before it; ff_idx unaffected (still earlier)
    out = out[:ff_idx + 1] + prm_block + out[ff_idx + 1:]
    out = out + mol_lines

    with open(args.out_topol, "w") as fh:
        fh.writelines(out)

    print(f"[OK] Wrote {args.out_topol}")
    print(f"     parameters inserted after line {ff_idx + 1}")
    print(f"     molecule topologies inserted before the water include")
    print(f"     [ molecules ] extended with: {', '.join(m['name'] for m in args.molecule)}")


if __name__ == "__main__":
    main()
