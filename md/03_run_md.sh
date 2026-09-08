#!/bin/bash --login
# System build + minimisation + NVT + NPT + production MD, after the
# protein and each molecule (ligand + GSH) are already prepped
# (01_prep_protein.sh, 02_prep_molecule.sh x2). Assembles topol.top,
# merges coordinates, solvates, neutralises (auto-detected SOL group),
# energy-minimises, builds restraints + the Protein_LIG_UNL index group,
# then NVT -> NPT -> MD.
#
# Validated end-to-end 2026-08-10 (GSTA1-1 + GSH + lig_03506): production
# confirmed at ~43 ns/day. Two real bugs fixed along the way -- see
# docs/implementation_notes.md.
#
# Expects (paths passed in, nothing hardcoded):
#   --protein-dir   output of 01_prep_protein.sh
#   --gsh-dir       output of 02_prep_molecule.sh for GSH (RESNAME=LIG --
#                   see docs/implementation_notes.md for why)
#   --ligand-dir    output of 02_prep_molecule.sh for this ligand (RESNAME=UNL)
#   --mdp-dir       directory with ions.mdp/minim.mdp/nvt.mdp/npt.mdp/md.mdp
#   --ffdir         charmm36-jul2022.ff directory
#   --outdir        where this run's output goes
#
# Usage:
#   ./03_run_md.sh --protein-dir system/protein --gsh-dir system/gsh \
#       --ligand-dir runs/lig00001/prep --mdp-dir ../md/mdp \
#       --ffdir charmm36-jul2022.ff --outdir runs/lig00001
set -e

TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tools"

while [ $# -gt 0 ]; do
    case "$1" in
        --protein-dir) PROTEIN_DIR=$2; shift 2 ;;
        --gsh-dir) GSH_DIR=$2; shift 2 ;;
        --ligand-dir) LIGAND_DIR=$2; shift 2 ;;
        --mdp-dir) MDP_DIR=$2; shift 2 ;;
        --ffdir) FFDIR=$2; shift 2 ;;
        --outdir) OUTDIR=$2; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

for v in PROTEIN_DIR GSH_DIR LIGAND_DIR MDP_DIR FFDIR OUTDIR; do
    if [ -z "${!v}" ]; then
        flag=$(echo "$v" | tr '[:upper:]_' '[:lower:]-')
        echo "ERROR: --${flag} is required"
        echo "Usage: $0 --protein-dir D --gsh-dir D --ligand-dir D --mdp-dir D --ffdir D --outdir D"
        exit 1
    fi
done

for f in "$PROTEIN_DIR/protein_processed.gro" "$PROTEIN_DIR/topol.top" \
         "$GSH_DIR/lig.itp" "$GSH_DIR/lig.prm" "$GSH_DIR/lig_ini.pdb" \
         "$LIGAND_DIR/unl.itp" "$LIGAND_DIR/unl.prm" "$LIGAND_DIR/unl_ini.pdb" \
         "$MDP_DIR/ions.mdp" "$MDP_DIR/minim.mdp" "$MDP_DIR/nvt.mdp" \
         "$MDP_DIR/npt.mdp" "$MDP_DIR/md.mdp"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: expected input not found: $f"
        exit 1
    fi
done

# Resolve every input to absolute before cd'ing into OUTDIR below.
PROTEIN_DIR=$(cd "$PROTEIN_DIR" && pwd)
GSH_DIR=$(cd "$GSH_DIR" && pwd)
LIGAND_DIR=$(cd "$LIGAND_DIR" && pwd)
MDP_DIR=$(cd "$MDP_DIR" && pwd)
FFDIR_ABS=$(cd "$(dirname "$FFDIR")" && pwd)/$(basename "$FFDIR")

mkdir -p "$OUTDIR"
cd "$OUTDIR"

ln -sfn "$FFDIR_ABS" "./$(basename "$FFDIR")"
FFDIR_NAME="./$(basename "$FFDIR")"

echo "======================================================================"
echo "STAGE 1: assemble topol.top and complex coordinates"
echo "======================================================================"

cp "$GSH_DIR/lig.itp" "$GSH_DIR/lig.prm" .
cp "$LIGAND_DIR/unl.itp" "$LIGAND_DIR/unl.prm" .

# pdb2gmx's own position-restraint file (posre*.itp) must be copied in --
# topol.top's #ifdef POSRES block expects it locally, not just in PROTEIN_DIR.
shopt -s nullglob
POSRE_FILES=("$PROTEIN_DIR"/posre*.itp)
shopt -u nullglob
if [ ${#POSRE_FILES[@]} -eq 0 ]; then
    echo "ERROR: no posre*.itp found in $PROTEIN_DIR -- topol.top's #ifdef POSRES"
    echo "block needs one. Re-run 01_prep_protein.sh, or check it actually completed."
    exit 1
fi
cp "${POSRE_FILES[@]}" .

echo "[1/3] gmx_mpi editconf: GSH/UNL pdb -> gro"
gmx_mpi editconf -f "$GSH_DIR/lig_ini.pdb" -o gsh.gro
gmx_mpi editconf -f "$LIGAND_DIR/unl_ini.pdb" -o unl.gro

echo "[2/3] Building topol.top (replaces 'manually add itp/prm to topol.top')"
# molecule name here MUST be "LIG" -- matches the [ moleculetype ] name
# cgenff_charmm2gmx.py wrote inside lig.itp (from the .str's RESI LIG line).
python3 "$TOOLS_DIR/build_topol.py" "$PROTEIN_DIR/topol.top" topol.top \
    --molecule "LIG:lig.itp:lig.prm:posre_LIG.itp" \
    --molecule "UNL:unl.itp:unl.prm:posre_unl.itp"
# posre_LIG.itp / posre_unl.itp don't exist yet -- fine, wrapped in #ifdef
# POSRES; stage 3 below (genrestr) produces them before EM/NVT/NPT need them.

echo "[3/3] Merging protein + GSH + UNL coordinates (order MUST match topol.top's [ molecules ])"
python3 "$TOOLS_DIR/merge_gro.py" complex.gro \
    "$PROTEIN_DIR/protein_processed.gro" gsh.gro unl.gro

echo ""
echo "======================================================================"
echo "STAGE 2: solvate + neutralise"
echo "======================================================================"
gmx_mpi editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0
gmx_mpi solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

echo "[genion] auto-detecting the SOL group index (no hardcoded 'select group 15')"
gmx_mpi grompp -f "$MDP_DIR/ions.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
NDX_LISTING=$(printf "q\n" | gmx_mpi make_ndx -f solv.gro -o solv_default.ndx 2>&1)
SOL_IDX=$(echo "$NDX_LISTING" | grep -E "^\s*[0-9]+\s+SOL\s" | head -1 | awk '{print $1}')
if [ -z "$SOL_IDX" ]; then
    echo "ERROR: could not auto-detect the SOL group -- inspect solv_default.ndx / the listing below:"
    echo "$NDX_LISTING"
    exit 1
fi
echo "  Detected SOL group index: $SOL_IDX"
printf "%s\n" "$SOL_IDX" | gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

echo ""
echo "======================================================================"
echo "STAGE 3: energy minimisation + restraints"
echo "======================================================================"
gmx_mpi grompp -f "$MDP_DIR/minim.mdp" -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
gmx_mpi mdrun -v -deffnm em

echo "[make_ndx] non-hydrogen atoms for UNL restraints"
UNL_NDX_LISTING=$(printf "0 & ! a H*\nq\n" | gmx_mpi make_ndx -f unl.gro -o index_unl.ndx 2>&1)
UNL_NOH_IDX=$(echo "$UNL_NDX_LISTING" | grep -E "^\s*[0-9]+\s+\S+\s*:" | awk '{print $1}' | sort -n | tail -1)

echo "[make_ndx] non-hydrogen atoms for GSH restraints"
GSH_NDX_LISTING=$(printf "0 & ! a H*\nq\n" | gmx_mpi make_ndx -f gsh.gro -o index_GSH.ndx 2>&1)
GSH_NOH_IDX=$(echo "$GSH_NDX_LISTING" | grep -E "^\s*[0-9]+\s+\S+\s*:" | awk '{print $1}' | sort -n | tail -1)

if [ -z "$UNL_NOH_IDX" ] || [ -z "$GSH_NOH_IDX" ]; then
    echo "ERROR: could not determine the non-hydrogen group index for UNL/GSH restraints"
    exit 1
fi
echo "  UNL non-H group index: $UNL_NOH_IDX | GSH non-H group index: $GSH_NOH_IDX"

printf "%s\n" "$UNL_NOH_IDX" | gmx_mpi genrestr -f unl.gro -n index_unl.ndx -o posre_unl.itp -fc 1000 1000 1000
printf "%s\n" "$GSH_NOH_IDX" | gmx_mpi genrestr -f gsh.gro -n index_GSH.ndx -o posre_LIG.itp -fc 1000 1000 1000

echo "[make_ndx] combined Protein_LIG_UNL group for tc-grps (nvt.mdp/npt.mdp/md.mdp expect this exact name)"
NDX_LISTING2=$(printf "q\n" | gmx_mpi make_ndx -f em.gro -o index.ndx 2>&1)
PROT_IDX=$(echo "$NDX_LISTING2" | grep -E "^\s*[0-9]+\s+Protein\s" | head -1 | awk '{print $1}')
# GSH's default group is named "LIG" (its RESI name), not "GSH".
GSH_IDX=$(echo "$NDX_LISTING2" | grep -E "^\s*[0-9]+\s+LIG\s" | head -1 | awk '{print $1}')
UNL_IDX=$(echo "$NDX_LISTING2" | grep -E "^\s*[0-9]+\s+UNL\s" | head -1 | awk '{print $1}')
if [ -z "$PROT_IDX" ] || [ -z "$GSH_IDX" ] || [ -z "$UNL_IDX" ]; then
    echo "ERROR: could not find Protein/LIG(GSH)/UNL default groups in em.gro's index listing:"
    echo "$NDX_LISTING2"
    exit 1
fi
echo "  Protein=$PROT_IDX  GSH=$GSH_IDX  UNL=$UNL_IDX -- combining into Protein_LIG_UNL"
MAX_IDX=$(echo "$NDX_LISTING2" | grep -E "^\s*[0-9]+\s+\S+\s*:" | awk '{print $1}' | sort -n | tail -1)
NEW_IDX=$(( MAX_IDX + 1 ))
printf "%s | %s | %s\nname %s Protein_LIG_UNL\nq\n" \
    "$PROT_IDX" "$GSH_IDX" "$UNL_IDX" "$NEW_IDX" \
    | gmx_mpi make_ndx -f em.gro -n index.ndx -o index.ndx

echo ""
echo "======================================================================"
echo "STAGE 4: NVT -> NPT -> production MD"
echo "======================================================================"
echo "[NVT]"
gmx_mpi grompp -f "$MDP_DIR/nvt.mdp" -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 1
gmx_mpi mdrun -v -deffnm nvt

echo "[NPT]"
gmx_mpi grompp -f "$MDP_DIR/npt.mdp" -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 1
gmx_mpi mdrun -v -deffnm npt

echo "[MD]"
gmx_mpi grompp -f "$MDP_DIR/md.mdp" -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_5.tpr -maxwarn 1
gmx_mpi mdrun -deffnm md_0_5

echo ""
echo "[OK] Production MD complete: $OUTDIR/md_0_5.xtc / .gro / .edr / .log"
echo "Next: 04_analyze_md.sh (QC/RMSD/warhead-distance) and 05_run_mmpbsa.sh (MM-PBSA rescoring)"
