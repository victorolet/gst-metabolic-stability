#!/bin/bash --login
# =============================================================================
# 02_prep_molecule.sh -- ligand OR GSH preprocessing (same procedure for both,
# per Harry's doc: "As per 'Ligand Preprocessing' but renaming UNL to GSH")
# =============================================================================
# Automates:
#   - docked-pose PDBQT -> mol2 (obabel), or accepts a mol2 directly (GSH
#     co-crystal molecules typically start as mol2/pdb already)
#   - sort_mol2_bonds.pl
#   - residue renaming (fix_mol2_resname.py -- robust to whatever placeholder
#     name the upstream tool used, not just a literal '*****')
#
# Then, because CGenFF itself has no license set up yet (webserver-only for
# now), this script PAUSES here: it tells you where to submit the fixed
# mol2 and where to save the resulting .str file, then exits 0. Run the
# exact same command again afterwards -- it detects the .str file is now
# present and finishes the job (run cgenff_charmm2gmx.py), producing
# <resname_lc>.itp/.prm/.top/_ini.pdb. See the LP/LPH note further down for
# why this does NOT strip LP/LPH lines by default, unlike Harry's original doc.
#
# Usage:
#   ./02_prep_molecule.sh input.pdbqt|input.mol2 RESNAME workdir/ [ffdir]
#   STRIP_LP=1 ./02_prep_molecule.sh ...   # force the old LP/LPH-stripping behaviour
#
# Examples:
#   ./02_prep_molecule.sh docked_pose_lig00001.pdbqt UNL runs/lig00001/
#   ./02_prep_molecule.sh gsh_cofactor.mol2 GSH system/gsh/
# =============================================================================
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
mkdir -p "$WORKDIR"
RESNAME_LC=$(echo "$RESNAME" | tr '[:upper:]' '[:lower:]')
STR_FILE="$WORKDIR/${RESNAME}.str"
FIXED_MOL2="$WORKDIR/${RESNAME}_fixed.mol2"

# ---- Stage 1: mol2 prep (cheap, always re-run so edits to INPUT propagate) ----
case "$INPUT" in
    *.pdbqt)
        echo "[1/3] Converting docked pose PDBQT -> mol2 (obabel)"
        obabel "$INPUT" -O "$WORKDIR/raw.mol2" 2>&1 | tail -5
        ;;
    *.mol2)
        cp "$INPUT" "$WORKDIR/raw.mol2"
        ;;
    *)
        echo "ERROR: input must be .pdbqt or .mol2, got: $INPUT"
        exit 1
        ;;
esac

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

# NOTE ON LP/LPH: Harry's original doc says to manually strip LP/LPH lines
# from the .str before conversion. This version of cgenff_charmm2gmx.py
# (fetched fresh from Lemkul-Lab, see tools/cgenff_charmm2gmx.py header) has
# built-in lone-pair support added specifically for CGenFF >=4.0 halogens
# (LONE parsing, is_lp(), 2fd virtual-site construction) -- stripping LP
# lines ahead of it would likely disable correct handling of halogenated
# warheads (exactly the SNAr/electrophile classes this project's alert_score
# descriptor flags), which the script is designed to handle properly on its
# own. So by default this script does NOT strip LP/LPH lines. If Harry
# confirms his version of the workflow still needs that (e.g. an older
# converter script without LP support), pass --strip-lp to force the old
# behaviour -- but verify against a known-good molecule first either way.
STR_TO_USE="$STR_FILE"
if [ "$STRIP_LP" = "1" ]; then
    echo "[--strip-lp set] Removing LP/LPH lines from the stream file"
    grep -vE '^\s*(ATOM\s+LP|LONE\s)' "$STR_FILE" > "$WORKDIR/${RESNAME}_clean.str"
    STR_TO_USE="$WORKDIR/${RESNAME}_clean.str"
fi

echo "Running cgenff_charmm2gmx.py"
cd "$WORKDIR"
python3 "$TOOLS_DIR/cgenff_charmm2gmx.py" "$RESNAME" \
    "$(basename "$FIXED_MOL2")" "$(basename "$STR_TO_USE")" "$FFDIR"

echo ""
echo "[OK] Molecule prep complete: ${RESNAME_LC}.itp / ${RESNAME_LC}.prm / ${RESNAME_LC}_ini.pdb in $WORKDIR"
