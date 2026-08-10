#!/usr/bin/env python3
# =============================================================================
# check_warhead_proximity.py -- automatable proxy for "does this docked pose
# make chemical sense", scaled across a whole batch of docking results.
# =============================================================================
#
# What this checks: for each docked ligand pose, finds the ligand's
# GSH-reactive electrophilic "warhead" atom (the same structural-alert
# definition the ML pipeline uses to flag reactive centres) and measures its
# 3D distance to GSH's thiolate sulfur (SG2) in the receptor. A pose whose
# warhead sits close to SG2 is plausible for the actual conjugation
# chemistry this project cares about; one where the warhead points away
# from SG2 -- even with a good Vina score -- probably isn't, since Vina's
# score reflects total binding energy, not reactive-geometry.
#
# This does NOT replace visual inspection entirely -- it can't catch every
# kind of implausible pose (clashes, strained geometry) the way eyeballing
# it in PyMOL/ChimeraX can. Treat it as a cheap triage filter: run it across
# a whole results/ batch to rank/flag candidates, then visually confirm
# whichever one(s) you're about to commit real MD compute to.
#
# ELECTRO_SMARTS below is copied from ml_pipeline/pipeline_HAL20260712_
# ensemble.py's own ELECTRO_SMARTS dict (not imported -- that module pulls
# in sklearn/matplotlib/imblearn/tqdm as import-time dependencies, heavy
# for what should be a lightweight geometry check). KEEP THESE IN SYNC
# MANUALLY if the ML pipeline's structural-alert list ever changes --
# that's the one real maintenance cost of not importing it directly.
#
# How atom correspondence works: meeko (which produced these ligand PDBQTs
# via prepare_ligands.py) writes "REMARK SMILES <smiles>" plus one or more
# "REMARK SMILES IDX <s1> <p1> <s2> <p2> ..." lines into every PDBQT it
# writes, pairing each 1-indexed SMILES heavy-atom index with its PDBQT
# atom serial number. This script re-parses THAT embedded SMILES (not the
# original training_data.csv SMILES -- RDKit/meeko may reorder atoms during
# processing, so re-deriving from the exact string meeko recorded is the
# only safe way to keep indices consistent), finds warhead atoms via
# ELECTRO_SMARTS, and uses the embedded index map to find the matching
# PDBQT atom's 3D coordinates. Verified against a real docked pose: the
# SMILES indices in the REMARK exactly span 1..N_heavy_atoms with no gaps,
# while PDBQT serials range higher (polar hydrogens get their own PDBQT
# atoms but aren't part of the heavy-atom-only SMILES numbering).
#
# Usage:
#   python3 check_warhead_proximity.py --receptor receptors/1PKW_GSH.pdbqt \
#       --results-glob "results/*/*.pdbqt" --out warhead_proximity.csv
#
# Requires: rdkit (already in the gst_ml env)
# =============================================================================

import argparse
import csv
import glob
import math
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# ---- copied from ml_pipeline/pipeline_HAL20260712_ensemble.py -- keep in sync ----
ELECTRO_SMARTS = {
    "snar_p_no2": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:c[N+](=O)[O-]"), 2.0),
    "snar_p_cn": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:cC#N"), 1.0),
    "snar_p_cf3": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:cC(F)(F)F"), 1.0),
    "snar_p_sulfone": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:cS(=O)(=O)[C,c,N]"), 1.0),
    "snar_p_carbonyl": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:cC(=O)"), 1.0),
    "snar_p_ammonium": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:c[N+]"), 1.0),
    "snar_p_noxide": (Chem.MolFromSmarts("[n+]1([O-])ccc([F,Cl,Br])cc1"), 1.0),
    "snar_o_no2": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c[N+](=O)[O-]"), 2.0),
    "snar_o_cn": (Chem.MolFromSmarts("[c]([F,Cl,Br]):cC#N"), 1.0),
    "snar_o_cf3": (Chem.MolFromSmarts("[c]([F,Cl,Br]):cC(F)(F)F"), 1.0),
    "snar_o_sulfone": (Chem.MolFromSmarts("[c]([F,Cl,Br]):cS(=O)(=O)[C,c,N]"), 1.0),
    "snar_o_carbonyl": (Chem.MolFromSmarts("[c]([F,Cl,Br]):cC(=O)"), 1.0),
    "snar_o_ammonium": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c[N+]"), 1.0),
    "snar_o_noxide": (Chem.MolFromSmarts("[n+]1([O-])c([F,Cl,Br])cccc1"), 1.0),
    "snar_2_pyridyl": (Chem.MolFromSmarts("[c]([F,Cl,Br]):n"), 1.5),
    "snar_3_pyridyl": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:n"), 1.5),
    "snar_pyridazine": (Chem.MolFromSmarts("[c]([F,Cl,Br])(:n):n"), 1.5),
    "snar_pyrimidine": (Chem.MolFromSmarts("[c]([F,Cl,Br]):n:c:n"), 1.5),
    "snar_trazine": (Chem.MolFromSmarts("[c]([F,Cl,Br]):n:c:n:c:n"), 1.5),
    "snar_o_trihalo": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c[F,Cl,Br,I]"), 1.0),
    "snar_p_trihalo": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:c[F,Cl,Br,I]"), 1.0),
    "snar_sulfonate_LG1": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:c[N+](=O)[O-]"), 1.0),
    "snar_sulfonate_LG2": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:cC#N"), 1.0),
    "snar_sulfonate_LG3": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:cC(F)(F)F"), 1.0),
    "snar_sulfonate_LG4": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:cC(=O)"), 1.0),
    "snar_sulfonate_LG5": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:c[N+]"), 1.0),
    "snar_sulfonate_LG6": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c[N+](=O)[O-]"), 1.0),
    "snar_sulfonate_LG7": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):cC#N"), 1.0),
    "snar_sulfonate_LG8": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):cC(=O)"), 1.0),
    "snar_sulfonate_LG9": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):n"), 1.0),
    "snar_sulfonate_LG10": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:n"), 1.0),
    "snar_sulfone_LG1": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):c:c:c[N+](=O)[O-]"), 1.0),
    "snar_sulfone_LG2": (Chem.MolFromSmarts("[c](S(=O)(=O)c):c:c:c[N+](=O)[O-]"), 1.0),
    "snar_sulfone_LG3": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):c:c:cC#N"), 1.0),
    "snar_sulfone_LG4": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):c:c:cC(=O)"), 1.0),
    "snar_sulfone_LG5": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):c[N+](=O)[O-]"), 1.0),
    "snar_sulfone_LG6": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):cC#N"), 1.0),
    "snar_sulfone_LG7": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):cC(=O)"), 1.0),
    "snar_sulfone_LG8": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):n"), 1.0),
    "snar_sulfone_LG9": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):c:c:n"), 1.0),
    "ab_unsat_carbonyl": (Chem.MolFromSmarts("[CX3:1]=[CX3]C=O"), 1.0),
    "acrylate_ester": (Chem.MolFromSmarts("[CX3:1]=[CX3]C(=O)[O,N]"), 1.0),
    "acrylonitrile": (Chem.MolFromSmarts("[CX3:1]=[CX3]C#N"), 1.0),
    "vinyl_sulfone": (Chem.MolFromSmarts("[CX3:1]=[CX3]S(=O)(=O)[C,c]"), 1.0),
    "nitroalkene": (Chem.MolFromSmarts("[CX3:1]=[CX3][N+](=O)[O-]"), 1.0),
    "divinyl_enone": (Chem.MolFromSmarts("[CX3:1]=[CX3]C(=O)C"), 1.0),
    "prim_ahalide": (Chem.MolFromSmarts("[CX4;!$(C(F)(F)F):1][Cl,Br,I]"), 1.0),
    "asulfonate": (Chem.MolFromSmarts("[CX4:1]OS(=O)(=O)[C,c]"), 1.0),
    "atriflate": (Chem.MolFromSmarts("[CX4:1]OS(=O)(=O)C(F)(F)F"), 1.0),
    "aammonium": (Chem.MolFromSmarts("[CH2X4:1][N+](C)(C)C"), 1.0),
    "epoxide": (Chem.MolFromSmarts("[CX4:1]1[OX2][CX4]1"), 1.5),
}
# ---- end copied block ----


def find_warhead_atom_indices(smiles):
    """0-indexed RDKit atom indices of every structural-alert match's first atom."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    idxs = set()
    for patt, _weight in ELECTRO_SMARTS.values():
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            if match:
                idxs.add(match[0])
    return sorted(idxs)


def parse_pdbqt_ligand(path):
    """Return (smiles, {rdkit_0idx: pdbqt_serial}, {pdbqt_serial: (x,y,z)})
    for the first MODEL (Vina's mode 1 = best pose) -- or the whole file if
    there's no MODEL/ENDMDL wrapper at all (single-pose PDBQTs)."""
    smiles = None
    idx_map = {}
    coords = {}
    model_seen = 0
    in_model_1 = True
    with open(path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                model_seen += 1
                in_model_1 = (model_seen == 1)
                continue
            if line.startswith("ENDMDL"):
                if in_model_1:
                    break
                continue
            if not in_model_1:
                continue
            if line.startswith("REMARK SMILES IDX"):
                nums = [int(x) for x in line.split()[3:]]
                for i in range(0, len(nums) - 1, 2):
                    smiles_idx_1based, pdbqt_serial = nums[i], nums[i + 1]
                    idx_map[smiles_idx_1based - 1] = pdbqt_serial
            elif line.startswith("REMARK SMILES") and not line.startswith("REMARK SMILES IDX"):
                smiles = line.split(None, 2)[2].strip()
            elif line.startswith(("ATOM", "HETATM")):
                serial = int(line[6:11])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords[serial] = (x, y, z)
    return smiles, idx_map, coords


def find_atom_coords(pdb_path, resname, atomname):
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                aname = line[12:16].strip()
                rname = line[17:20].strip()
                if rname == resname and aname == atomname:
                    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return None


def dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--receptor", required=True, help="Receptor PDBQT with GSH bound (e.g. receptors/1PKW_GSH.pdbqt)")
    ap.add_argument("--results-glob", required=True, help='Glob for docked ligand PDBQTs, e.g. "results/*/*.pdbqt"')
    ap.add_argument("--out", default="warhead_proximity.csv")
    ap.add_argument("--sg-resname", default="GSH")
    ap.add_argument("--sg-atomname", default="SG2", help="GSH's thiolate sulfur atom name (see receptor HETATM records)")
    args = ap.parse_args()

    sg2 = find_atom_coords(args.receptor, args.sg_resname, args.sg_atomname)
    if sg2 is None:
        sys.exit(f"ERROR: could not find {args.sg_resname} {args.sg_atomname} in {args.receptor}")

    paths = sorted(glob.glob(args.results_glob))
    if not paths:
        sys.exit(f"ERROR: no files matched --results-glob '{args.results_glob}'")

    rows = []
    for path in paths:
        smiles, idx_map, coords = parse_pdbqt_ligand(path)
        if not smiles:
            rows.append({"file": path, "distance_to_SG2": "", "status": "no REMARK SMILES found"})
            continue
        warhead_idxs = find_warhead_atom_indices(smiles)
        if not warhead_idxs:
            rows.append({"file": path, "distance_to_SG2": "", "status": "no structural-alert warhead matched"})
            continue
        best = None
        for wi in warhead_idxs:
            serial = idx_map.get(wi)
            if serial is None or serial not in coords:
                continue
            d = dist(coords[serial], sg2)
            if best is None or d < best:
                best = d
        if best is None:
            rows.append({"file": path, "distance_to_SG2": "", "status": "warhead atom not resolved to a PDBQT atom"})
        else:
            rows.append({"file": path, "distance_to_SG2": round(best, 2), "status": "ok"})

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["file", "distance_to_SG2", "status"])
        w.writeheader()
        w.writerows(rows)

    ok = sorted((r for r in rows if r["status"] == "ok"), key=lambda r: r["distance_to_SG2"])
    print(f"Wrote {args.out}: {len(rows)} poses checked, {len(ok)} with a resolved warhead-to-SG2 distance")
    print("(rough rule of thumb -- not a hard cutoff: a few A suggests a plausible pre-reactive")
    print(" orientation; tens of A means the warhead is pointed well away from SG2 regardless")
    print(" of Vina score. Always eyeball whichever pose you actually commit MD compute to.)")
    print()
    print("Closest 10 warhead-to-SG2 distances:")
    for r in ok[:10]:
        print(f"  {r['distance_to_SG2']:>6} A  {r['file']}")


if __name__ == "__main__":
    main()
