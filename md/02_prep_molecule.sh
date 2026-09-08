#!/bin/bash --login
# Ligand OR GSH prep (same procedure for both): mol2 conversion, bond
# sorting, residue rename, pause for manual CGenFF webserver submission,
# then cgenff_charmm2gmx.py -> itp/prm. RESNAME must match the .str file's
# internal RESI line exactly -- GSH's stream files say "RESI LIG", not
# "GSH"; docked-ligand files say "RESI UNL". Check with `grep '^RESI'
# file.str` first. Rationale + gotchas: docs/implementation_notes.md.
#
# Usage:
#   ./02_prep_molecule.sh input.pdbqt|input.pdb|input.mol2 RESNAME workdir/ [ffdir]
#   STRIP_LP=1 ./02_prep_molecule.sh ...   # force the old LP/LPH-stripping behaviour
#
# Examples:
#   ./02_prep_molecule.sh docked_pose_lig00001.pdbqt UNL runs/lig00001/
#   ./02_prep_molecule.sh gsh_from_1pkw.pdb LIG system/gsta1_gsh_prepped/
set -e

INPUT=$1
RESNAME=$2
WORKDIR=$3
FFDIR=${4:-charmm36-jul2022.ff}
STRIP_LP=${STRIP_LP:-0}

if [ -z "$INPUT" ] || [ -z "$RESNAME" ] || [ -z "$WORKDIR" ]; then
    echo "Usage: $0 input.pdbqt|input.mol2 RESNAME workdir/ [ffdir]"
    exit 1
fi
if [ ! -f "$INPUT" ]; then
    echo "ERROR: $INPUT not found"
    exit 1
fi

TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tools"
# Resolve FFDIR to absolute now, before cd'ing into WORKDIR below.
FFDIR_ABS=$(cd "$(dirname "$FFDIR")" && pwd)/$(basename "$FFDIR")
mkdir -p "$WORKDIR"
RESNAME_LC=$(echo "$RESNAME" | tr '[:upper:]' '[:lower:]')
STR_FILE="$WORKDIR/${RESNAME}.str"
FIXED_MOL2="$WORKDIR/${RESNAME}_fixed.mol2"

# ---- Stage 1: mol2 prep (cheap, always re-run so edits to INPUT propagate) ----
resolve_obabel() {
    # obabel isn't in either conda env (spack build) -- override > PATH > known path.
    SPACK_OBABEL="/software/projects/pawsey1376/volet/setonix/2025.08/software/linux-sles15-zen3/gcc-14.2.0/openbabel-3.1.1-onem6ip5yyp7artla3orqgwac5biq3ck/bin/obabel"
    if [ -n "$OBABEL_BIN" ]; then
        :
    elif command -v obabel >/dev/null 2>&1; then
        OBABEL_BIN=obabel
    elif [ -x "$SPACK_OBABEL" ]; then
        OBABEL_BIN="$SPACK_OBABEL"
    else
        echo "ERROR: obabel not found on PATH, no OBABEL_BIN override set, and"
        echo "the expected spack install isn't at: $SPACK_OBABEL"
        echo "Set OBABEL_BIN=/path/to/obabel and re-run."
        exit 1
    fi
}

case "$INPUT" in
    *.pdbqt)
        resolve_obabel
        echo "[1/3] Converting docked pose PDBQT -> mol2 ($OBABEL_BIN)"
        "$OBABEL_BIN" "$INPUT" -O "$WORKDIR/raw.mol2" 2>&1 | tail -5
        ;;
    *.pdb)
        resolve_obabel
        echo "[1/3] Converting PDB -> mol2 ($OBABEL_BIN)"
        "$OBABEL_BIN" "$INPUT" -O "$WORKDIR/raw.mol2" 2>&1 | tail -5
        ;;
    *.mol2)
        cp "$INPUT" "$WORKDIR/raw.mol2"
        ;;
    *)
        echo "ERROR: input must be .pdbqt, .pdb, or .mol2, got: $INPUT"
        exit 1
        ;;
esac

# Strip anything before the first TRIPOS line -- PyMOL-exported mol2s carry
# a leading comment line that sort_mol2_bonds.pl chokes on.
awk '/TRIPOS/{f=1} f' "$WORKDIR/raw.mol2" > "$WORKDIR/raw.mol2.tmp" && mv "$WORKDIR/raw.mol2.tmp" "$WORKDIR/raw.mol2"

echo "[2/3] Sorting mol2 bond order (sort_mol2_bonds.pl)"
perl "$TOOLS_DIR/sort_mol2_bonds.pl" "$WORKDIR/raw.mol2" "$WORKDIR/sorted.mol2"

echo "[3/3] Renaming residue -> $RESNAME (fix_mol2_resname.py)"
python3 "$TOOLS_DIR/fix_mol2_resname.py" "$WORKDIR/sorted.mol2" "$FIXED_MOL2" "$RESNAME"

# ---- Stage 2: CGenFF hand-off / resume ----
if [ ! -f "$STR_FILE" ]; then
    echo ""
    echo "=============================================================="
    echo "PAUSED -- CGenFF step needs to happen manually (no local CGenFF"
    echo "license set up yet; this is the one step not yet automated)."
    echo ""
    echo "  1. Submit this file to https://cgenff.com :"
    echo "       $FIXED_MOL2"
    echo "  2. Save the resulting stream file as:"
    echo "       $STR_FILE"
    echo "  3. Re-run this exact command -- it will detect the .str file"
    echo "     and finish automatically (run cgenff_charmm2gmx.py)."
    echo "=============================================================="
    exit 0
fi

echo ""
echo "Found $STR_FILE -- finishing molecule prep"

# LP/LPH lines are NOT stripped by default -- this cgenff_charmm2gmx.py
# build has native lone-pair support for CGenFF>=4.0 halogens; stripping
# would disable it. Pass --strip-lp for the old manual-doc behaviour.
STR_TO_USE="$STR_FILE"
if [ "$STRIP_LP" = "1" ]; then
    echo "[--strip-lp set] Removing LP/LPH lines from the stream file"
    grep -vE '^\s*(ATOM\s+LP|LONE\s)' "$STR_FILE" > "$WORKDIR/${RESNAME}_clean.str"
    STR_TO_USE="$WORKDIR/${RESNAME}_clean.str"
fi

echo "Running cgenff_charmm2gmx.py"
cd "$WORKDIR"
python3 "$TOOLS_DIR/cgenff_charmm2gmx.py" "$RESNAME" \
    "$(basename "$FIXED_MOL2")" "$(basename "$STR_TO_USE")" "$FFDIR_ABS"

echo ""
echo "[OK] Molecule prep complete: ${RESNAME_LC}.itp / ${RESNAME_LC}.prm / ${RESNAME_LC}_ini.pdb in $WORKDIR"
