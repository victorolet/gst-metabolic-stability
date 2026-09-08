#!/usr/bin/env python3
"""Identify the electrophilic 'warhead' atom in a ligand mol2 file (the same
CGenFF-ready mol2 used to build its GROMACS itp/prm), and print its 1-based
LOCAL atom number (position within this mol2's own @<TRIPOS>ATOM block) to
stdout.

Why a NUMBER and not the atom NAME: checked against a real ligand mol2 from
this project (md/ligand_mol2_cache/GSTA/3506_A_fixed.mol2) and found the
atom names are bare, non-unique element symbols ("C", "N", "H" repeated
25/5/27 times) -- unlike GSH's mol2, which does have unique per-atom names
(SG2, OE1, N1, ...). A name-based make_ndx selection ("UNL & a C") would
therefore select every carbon in the ligand, not the one warhead atom. A
1-based local atom number is unambiguous instead: 04_analyze_md.sh adds the
preceding protein+GSH atom counts to convert it into a whole-system atom
number valid in make_ndx's "a <N>" selector (merge_gro.py's fixed ordering
is protein, then GSH, then ligand, so this offset is deterministic).

ELECTRO_SMARTS below is copied from docking/check_warhead_proximity.py
(itself copied from ml_pipeline/pipeline_HAL20260712_ensemble.py) -- KEEP
THESE IN SYNC MANUALLY if the structural-alert list ever changes.

Why work from the mol2 directly instead of reusing check_warhead_proximity.py's
PDBQT/REMARK-SMILES-IDX approach: the mol2 used for CGenFF/MD prep is a
separately-prepared file from the docking PDBQT (different atom numbering),
so there's no safe index correspondence between the two. Matching SMARTS
directly against the mol2's own RDKit Mol keeps everything in one file's own
atom order (RDKit preserves mol2 atom order -- it doesn't renumber/
canonicalise), avoiding any cross-file index mapping.

Verified against a real project file (2026-08-13): correctly flagged atom
#56 of 3506_A_fixed.mol2 (a C.2 carbon, alpha to a C=C double bond feeding
a carbonyl -- an ab_unsat_carbonyl / Michael-acceptor-type match), a
chemically sensible result.

Usage: find_warhead_atom.py ligand_fixed.mol2
Prints: the winning atom's 1-based local atom number (e.g. "56") to stdout,
nothing else -- diagnostics (index, element, mol2 name) go to stderr. Exits
non-zero on any failure.
"""
import sys
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# ---- copied from docking/check_warhead_proximity.py -- keep in sync ----
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


def parse_mol2_atom_names(path):
    """Atom names in @<TRIPOS>ATOM file order (0-indexed), matching RDKit's
    own atom order -- RDKit's mol2 parser preserves atom order, it doesn't
    canonicalise/renumber. Used here only for the diagnostic message; the
    atom NUMBER (not the name) is what actually gets returned to the caller,
    since these names are frequently non-unique (see module docstring)."""
    names = []
    in_atoms = False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if s.startswith("@<TRIPOS>"):
                in_atoms = False
                continue
            if in_atoms and s:
                names.append(s.split()[1])
    return names


def find_warhead_atom_idx(mol):
    """0-indexed RDKit atom idx of the highest-weight structural-alert match's
    first atom, or None if nothing matched."""
    best = None  # (weight, atom_idx)
    for patt, weight in ELECTRO_SMARTS.values():
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            if match and (best is None or weight > best[0]):
                best = (weight, match[0])
    return best[1] if best else None


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: find_warhead_atom.py ligand_fixed.mol2")
    path = sys.argv[1]

    mol = Chem.MolFromMol2File(path, removeHs=False)
    if mol is None:
        sys.exit(f"ERROR: RDKit could not parse {path} (check bond typing / valences)")

    idx = find_warhead_atom_idx(mol)
    if idx is None:
        sys.exit("ERROR: no electrophilic/structural-alert warhead matched this ligand -- "
                 "check ELECTRO_SMARTS coverage, or this molecule may genuinely have none")

    names = parse_mol2_atom_names(path)
    if idx >= len(names):
        sys.exit(f"ERROR: RDKit atom index {idx} out of range for {len(names)} mol2 atom "
                 f"names -- RDKit/mol2 atom count or order diverged, inspect manually")

    local_atom_num = idx + 1  # 1-based, matches the mol2/unl.gro atom position
    element = mol.GetAtomWithIdx(idx).GetSymbol()
    print(f"[find_warhead_atom] {path}: RDKit atom idx {idx} (0-based) -> "
          f"local atom #{local_atom_num} ({element}, mol2 name '{names[idx]}')", file=sys.stderr)
    print(local_atom_num)


if __name__ == "__main__":
    main()
