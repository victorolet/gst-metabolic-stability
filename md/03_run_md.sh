#!/bin/bash --login
# =============================================================================
# 03_run_md.sh -- system build + minimisation + NVT + NPT + production MD
# =============================================================================
# Runs everything from Harry's doc AFTER the protein and each molecule
# (ligand + GSH) have already been prepped (01_prep_protein.sh,
# 02_prep_molecule.sh x2). Automates: topol.top assembly, complex
# construction, solvation, ion neutralisation (with automatic SOL group
# detection -- no hardcoded "select group 15" like the manual workflow),
# energy minimisation, restraint generation, the combined
# Protein_LIG_UNL index group the .mdp files expect, then NVT -> NPT -> MD.
#
# This has NOT been run end-to-end against a real protein/ligand system --
# there's no test data available yet to validate it against. Treat the
# first real run as a validation run: watch each stage's log for the usual
# GROMACS sanity signs (LINCS warnings, temperature/pressure drift, etc.)
# before trusting the output.
#
# Expects (paths passed in, nothing hardcoded):
#   --protein-dir   output of 01_prep_protein.sh (protein_processed.gro, topol.top)
#   --gsh-dir       output of 02_prep_molecule.sh for GSH (gsh.itp/.prm/_ini.pdb)
#   --ligand-dir    output of 02_prep_molecule.sh for this ligand (unl.itp/.prm/_ini.pdb)
#   --mdp-dir       directory with ions.mdp/minim.mdp/nvt.mdp/npt.mdp/md.mdp
#   --ffdir         charmm36-jul2022.ff directory
#   --outdir        where this run's output goes (e.g. runs/lig00001/)
#
# Usage:
#   ./03_run_md.sh --protein-dir system/protein --gsh-dir system/gsh \
#       --ligand-dir runs/lig00001/prep --mdp-dir ../md/mdp \
#       --ffdir charmm36-jul2022.ff --outdir runs/lig00001
# =============================================================================
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
         "$GSH_DIR/gsh.itp" "$GSH_DIR/gsh.prm" "$GSH_DIR/gsh_ini.pdb" \
         "$LIGAND_DIR/unl.itp" "$LIGAND_DIR/unl.prm" "$LIGAND_DIR/unl_ini.pdb" \
         "$MDP_DIR/ions.mdp" "$MDP_DIR/minim.mdp" "$MDP_DIR/nvt.mdp" \
         "$MDP_DIR/npt.mdp" "$MDP_DIR/md.mdp"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: expected input not found: $f"
        exit 1
    fi
done

mkdir -p "$OUTDIR"
cd "$OUTDIR"

FFDIR_ABS=$(cd "$(dirname "$FFDIR")" && pwd)/$(basename "$FFDIR")
ln -sfn "$FFDIR_ABS" "./$(basename "$FFDIR")"
FFDIR_NAME="./$(basename "$FFDIR")"

echo "======================================================================"
echo "STAGE 1: assemble topol.top and complex coordinates"
echo "======================================================================"

cp "$GSH_DIR/gsh.itp" "$GSH_DIR/gsh.prm" .
cp "$LIGAND_DIR/unl.itp" "$LIGAND_DIR/unl.prm" .

echo "[1/3] gmx editconf: GSH/UNL pdb -> gro"
gmx editconf -f "$GSH_DIR/gsh_ini.pdb" -o gsh.gro
gmx editconf -f "$LIGAND_DIR/unl_ini.pdb" -o unl.gro

echo "[2/3] Building topol.top (replaces 'manually add itp/prm to topol.top')"
python3 "$TOOLS_DIR/build_topol.py" "$PROTEIN_DIR/topol.top" topol.top \
    --molecule "GSH:gsh.itp:gsh.prm:posre_GSH.itp" \
    --molecule "UNL:unl.itp:unl.prm:posre_unl.itp"
# NOTE: posre_GSH.itp / posre_unl.itp don't exist yet at this point -- that's
# fine, they're wrapped in #ifdef POSRES and only EM/NVT/NPT (which do
# define POSRES) will need them, by which point genrestr (stage 3 below)
# has produced them.

echo "[3/3] Merging protein + GSH + UNL coordinates (order MUST match topol.top's [ molecules ])"
python3 "$TOOLS_DIR/merge_gro.py" complex.gro \
    "$PROTEIN_DIR/protein_processed.gro" gsh.gro unl.gro

echo ""
echo "======================================================================"
echo "STAGE 2: solvate + neutralise"
echo "======================================================================"
gmx editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

echo "[genion] auto-detecting the SOL group index (no hardcoded 'select group 15')"
gmx grompp -f "$MDP_DIR/ions.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
NDX_LISTING=$(printf "q\n" | gmx make_ndx -f solv.gro -o solv_default.ndx 2>&1)
SOL_IDX=$(echo "$NDX_LISTING" | grep -E "^\s*[0-9]+\s+SOL\s" | head -1 | awk '{print $1}')
if [ -z "$SOL_IDX" ]; then
    echo "ERROR: could not auto-detect the SOL group -- inspect solv_default.ndx / the listing below:"
    echo "$NDX_LISTING"
    exit 1
fi
echo "  Detected SOL group index: $SOL_IDX"
printf "%s\n" "$SOL_IDX" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

echo ""
echo "======================================================================"
echo "STAGE 3: energy minimisation + restraints"
echo "======================================================================"
gmx grompp -f "$MDP_DIR/minim.mdp" -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em

echo "[make_ndx] non-hydrogen atoms for UNL restraints"
UNL_NDX_LISTING=$(printf "0 & ! a H*\nq\n" | gmx make_ndx -f unl.gro -o index_unl.ndx 2>&1)
# the custom group we just made is always the highest-numbered group in the
# listing (make_ndx appends new groups after every pre-existing default one)
UNL_NOH_IDX=$(echo "$UNL_NDX_LISTING" | grep -E "^\s*[0-9]+\s+\S+\s*:" | awk '{print $1}' | sort -n | tail -1)

echo "[make_ndx] non-hydrogen atoms for GSH restraints"
GSH_NDX_LISTING=$(printf "0 & ! a H*\nq\n" | gmx make_ndx -f gsh.gro -o index_GSH.ndx 2>&1)
GSH_NOH_IDX=$(echo "$GSH_NDX_LISTING" | grep -E "^\s*[0-9]+\s+\S+\s*:" | awk '{print $1}' | sort -n | tail -1)

if [ -z "$UNL_NOH_IDX" ] || [ -z "$GSH_NOH_IDX" ]; then
    echo "ERROR: could not determine the non-hydrogen group index for UNL/GSH restraints"
    exit 1
fi
echo "  UNL non-H group index: $UNL_NOH_IDX | GSH non-H group index: $GSH_NOH_IDX"

printf "%s\n" "$UNL_NOH_IDX" | gmx genrestr -f unl.gro -n index_unl.ndx -o posre_unl.itp -fc 1000 1000 1000
printf "%s\n" "$GSH_NOH_IDX" | gmx genrestr -f gsh.gro -n index_GSH.ndx -o posre_GSH.itp -fc 1000 1000 1000
# (posre_unl.itp / posre_GSH.itp now exist -- topol.top's #ifdef POSRES
# blocks referencing them, written in stage 1, are now satisfiable.)

echo "[make_ndx] combined Protein_LIG_UNL group for tc-grps (nvt.mdp/npt.mdp/md.mdp expect this exact name)"
NDX_LISTING2=$(printf "q\n" | gmx make_ndx -f em.gro -o index.ndx 2>&1)
PROT_IDX=$(echo "$NDX_LISTING2" | grep -E "^\s*[0-9]+\s+Protein\s" | head -1 | awk '{print $1}')
GSH_IDX=$(echo "$NDX_LISTING2" | grep -E "^\s*[0-9]+\s+GSH\s" | head -1 | awk '{print $1}')
UNL_IDX=$(echo "$NDX_LISTING2" | grep -E "^\s*[0-9]+\s+UNL\s" | head -1 | awk '{print $1}')
if [ -z "$PROT_IDX" ] || [ -z "$GSH_IDX" ] || [ -z "$UNL_IDX" ]; then
    echo "ERROR: could not find Protein/GSH/UNL default groups in em.gro's index listing:"
    echo "$NDX_LISTING2"
    exit 1
fi
echo "  Protein=$PROT_IDX  GSH=$GSH_IDX  UNL=$UNL_IDX -- combining into Protein_LIG_UNL"
# the combined group make_ndx is about to create gets the next index after
# every group currently listed (max existing index + 1), not a guessed count
MAX_IDX=$(echo "$NDX_LISTING2" | grep -E "^\s*[0-9]+\s+\S+\s*:" | awk '{print $1}' | sort -n | tail -1)
NEW_IDX=$(( MAX_IDX + 1 ))
printf "%s | %s | %s\nname %s Protein_LIG_UNL\nq\n" \
    "$PROT_IDX" "$GSH_IDX" "$UNL_IDX" "$NEW_IDX" \
    | gmx make_ndx -f em.gro -n index.ndx -o index.ndx

echo ""
echo "======================================================================"
echo "STAGE 4: NVT -> NPT -> production MD"
echo "======================================================================"
echo "[NVT]"
gmx grompp -f "$MDP_DIR/nvt.mdp" -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 1
gmx mdrun -v -deffnm nvt

echo "[NPT]"
gmx grompp -f "$MDP_DIR/npt.mdp" -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 1
gmx mdrun -v -deffnm npt

echo "[MD]"
gmx grompp -f "$MDP_DIR/md.mdp" -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_5.tpr -maxwarn 1
gmx mdrun -deffnm md_0_5

echo ""
echo "[OK] Production MD complete: $OUTDIR/md_0_5.xtc / .gro / .edr / .log"
echo "Next: ligand RMSD (gmx rms, fit on protein backbone) and MM-GBSA rescoring"
echo "aren't in the doc Harry sent yet -- flag if you want those scripted too"
echo "once the exact reference/fit conventions are confirmed."
